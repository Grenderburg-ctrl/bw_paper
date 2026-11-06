# -*- coding: utf-8 -*-
"""
Парсер Государственного реестра медицинских изделий (Росздравнадзор).

Задача: найти ДЕЙСТВУЮЩИЕ регистрационные удостоверения (РУ) на наборы
для определения концентрации ИЛ-6, ИЛ-8, ИЛ-10 и NT-proBNP, пригодные
для исследования крови (цельная кровь / сыворотка / плазма).

Вместо Selenium + MHTML используется тот же JSON API, который вызывает
сам виджет https://elk.roszdravnadzor.gov.ru/widget/... :
  POST /public-gateway/registered-med-product/api/v1/med-product/filter-public   (поиск, все страницы)
  GET  /public-gateway/registered-med-product/api/v1/med-product/read-public/{id} (карточка РУ)
  GET  /public-gateway/registered-med-product/api/v1/model/read-by-product-id-public/{id} (модели/состав)

Сайт использует сертификат «Russian Trusted Root CA» (Минцифры), которого
нет в стандартном наборе certifi, поэтому скрипт скачивает официальные
корневой и промежуточный сертификаты с gu-st.ru и добавляет их к certifi.
Проверка TLS НЕ отключается.

Результат: pandas DataFrame `df_all` (все найденные действующие РУ с
разметкой) и `df` (только подходящие), плюс файл .xlsx.

В Google Colab: запустите ячейки по порядку (файл можно открыть как .ipynb
или вставить целиком в одну ячейку).
"""

# %% 0. Установка зависимостей (в Colab requests/pandas уже есть)
# !pip -q install openpyxl

# %% 1. Настройки
import os
import re
import time
import json
import certifi
import requests
import pandas as pd
from pathlib import Path

BASE = "https://elk.roszdravnadzor.gov.ru"
API = BASE + "/public-gateway/registered-med-product/api/v1"
ACTIVE_STATUS_ID = 1          # 1 = «Действует»
PAGE_SIZE = 100
DELAY = 0.3                   # пауза между запросами, сек
OUTPUT_XLSX = "RZN_IL6_IL8_IL10_NTproBNP.xlsx"
RAW_JSON = "rzn_raw_records.json"   # сырые ответы API (кэш)
USE_CACHE = False                   # True — не скачивать заново, взять RAW_JSON
# Чтобы сохранить на Google Drive:
# from google.colab import drive; drive.mount('/content/drive')
# OUTPUT_XLSX = "/content/drive/MyDrive/Colab Notebooks/RZN_IL6_IL8_IL10_NTproBNP.xlsx"

# Поисковые фразы (поиск сервера идёт по названию РУ + составу/моделям).
# Берётся объединение результатов, затем всё размечается локально.
SEARCH_TERMS = [
    "интерлейкин", "Interleukin", "цитокин",
    "IL-6", "IL6", "IL 6", "ИЛ-6", "ИЛ6",
    "IL-8", "IL8", "IL 8", "ИЛ-8", "ИЛ8",
    "IL-10", "IL10", "IL 10", "ИЛ-10", "ИЛ10",
    "NT-proBNP", "proBNP", "NT-pro", "про-BNP", "BNP",
    "натрийуретическ", "натрийуретического пептида",
]

# %% 2. HTTP-сессия с российским корневым сертификатом
CA_URLS = [
    "https://gu-st.ru/content/lending/russian_trusted_root_ca_pem.crt",
    "https://gu-st.ru/content/lending/russian_trusted_sub_ca_pem.crt",
]
CA_BUNDLE = Path("ca_bundle_with_russian_trusted_ca.pem")


def build_ca_bundle(path=CA_BUNDLE):
    parts = [Path(certifi.where()).read_text()]
    for url in CA_URLS:
        pem = requests.get(url, timeout=60).text.replace("\r", "").strip()
        if "BEGIN CERTIFICATE" not in pem:
            raise RuntimeError(f"Не удалось скачать сертификат: {url}")
        parts.append(pem)
    # каждый сертификат с новой строки (иначе OpenSSL выдаёт «PEM lib» ошибку)
    path.write_text("\n".join(p.strip() for p in parts) + "\n")
    return str(path)


CA_PATH = build_ca_bundle()
session = requests.Session()
session.headers.update({
    "User-Agent": "Mozilla/5.0 (research script)",
    "Accept": "application/json",
    "Referer": BASE + "/widget/",
})


def call(method, url, retries=8, **kw):
    """Запрос с повторами: сервер периодически рвёт соединения."""
    last = None
    for attempt in range(retries):
        try:
            # verify передаётся явно: session.verify перекрывается переменной REQUESTS_CA_BUNDLE
            r = session.request(method, url, timeout=90, verify=CA_PATH, **kw)
            if r.status_code == 404:
                return None
            r.raise_for_status()
            time.sleep(DELAY)
            return r.json()
        except requests.exceptions.SSLError as e:
            if "CERTIFICATE_VERIFY_FAILED" in str(e):
                raise
            last = e
            time.sleep(min(2 ** attempt, 30))
        except (requests.RequestException, ValueError) as e:
            last = e
            time.sleep(min(2 ** attempt, 30))
    raise RuntimeError(f"{method} {url} не удался после {retries} попыток: {last}")


# %% 3. Поиск по реестру (все страницы)
def search_all(term, status_ids=(ACTIVE_STATUS_ID,)):
    body = {"textSearch": term, "statusIds": list(status_ids)}
    rows, page = [], 0
    while True:
        data = call("POST", f"{API}/med-product/filter-public",
                    json=body, params={"page": page, "size": PAGE_SIZE})
        if not data:
            break
        rows.extend(data.get("content") or [])
        if data.get("last", True) or page + 1 >= data.get("totalPages", 0):
            break
        page += 1
    return rows


# %% 4. Карточки РУ и модели (состав набора)
def get_models(pid):
    names, page = [], 0
    while True:
        data = call("GET", f"{API}/model/read-by-product-id-public/{pid}",
                    params={"page": page, "size": 100, "sort": "id,asc"})
        if not data:
            break
        names += [(m.get("name") or "").strip() for m in data.get("content") or []]
        if data.get("last", True):
            break
        page += 1
    return names


def download_all():
    found = {}
    for term in SEARCH_TERMS:
        hits = search_all(term)
        for h in hits:
            found.setdefault(h["id"], {"hit": h, "terms": []})
            if term not in found[h["id"]]["terms"]:
                found[h["id"]]["terms"].append(term)
        print(f"{term!r:32} -> {len(hits):4d} действующих РУ")
    print(f"\nУникальных РУ: {len(found)}")

    details = {}
    for i, pid in enumerate(sorted(found), 1):
        card = call("GET", f"{API}/med-product/read-public/{pid}") or found[pid]["hit"]
        details[pid] = {"card": card, "models": get_models(pid), "terms": found[pid]["terms"]}
        if i % 20 == 0 or i == len(found):
            print(f"  загружено {i}/{len(found)}")
    return details


if USE_CACHE and Path(RAW_JSON).exists():
    details = {int(k): v for k, v in json.loads(Path(RAW_JSON).read_text(encoding="utf-8")).items()}
    print(f"Загружено из кэша {RAW_JSON}: {len(details)} РУ")
else:
    details = download_all()
    Path(RAW_JSON).write_text(json.dumps(details, ensure_ascii=False), encoding="utf-8")
    print("Сырые данные сохранены:", RAW_JSON)


# %% 5. Разметка: маркеры, тип образца, тип изделия
DASH = r"[\s\-‐-―−]*"
IL_PREFIX = r"(?:интерлейкин\w*|interleukin\w*|\bIL|\bИЛ|\bhIL|\bhuIL)"
# «IL-6», «интерлейкина 6», «IL-1β, -6, -8 и -10» и т.п.
IL_NUM = r"\d{1,2}(?![\d.,]\d)(?!\d)\s*[αβab]?"
IL_LIST = re.compile(
    IL_PREFIX + DASH + r"(" + IL_NUM +
    r"(?:\s*(?:,|/|;|\bи\b|\band\b|\bor\b|\bили\b)\s*-?\s*(?:IL" + DASH + r")?" + IL_NUM + r")*)",
    re.I)
# Рецепторы IL (sIL-6R, IL2R, «рецептор интерлейкина 6») — не считаем самим IL
IL_RECEPTOR = re.compile(r"рецептор\w*\s+(?:к\s+)?интерлейкин\w*" + DASH + r"\d+\w*"
                         r"|\bs?IL" + DASH + r"\d+\s*s?R(?:alpha|α|a)?\b", re.I)
NTPROBNP = re.compile(
    r"NT" + DASH + r"pro" + DASH + r"BNP|N" + DASH + r"терминальн\w*.{0,60}?натрийуретическ"
    r"|pro" + DASH + r"BNP|про" + DASH + r"BNP|предшественник\w*.{0,40}?натрийуретическ"
    r"|N" + DASH + r"terminal.{0,30}?pro.{0,10}?B" + DASH + r"type natriuretic", re.I | re.S)
BNP_ONLY = re.compile(r"(?<![\w-])BNP\b|мозгов\w*\s+натрийуретическ|натрийуретическ\w*\s+пептид\w*\s+(?:B|В)\b", re.I)

BLOOD = re.compile(r"кров|сыворот|плазм|serum|plasma|blood|капилляр", re.I)
OTHER_MATRIX = re.compile(r"моч[аеиу]|urine|ликвор|спинномозг|cerebrospinal|супернатант|культур\w* клет|"
                          r"cell culture|слюн|saliva|бронхоальвеол|экссудат|вагинальн|цервикальн|мазк|мазок|swab|"
                          r"амниотич|околоплодн|кал[ае]\b|feces|faeces", re.I)
# Вид-аналит (не «бычий альбумин» или «мышиные антитела» в составе набора)
ANIMAL = re.compile(
    r"(?:интерлейкин\w*|\bIL\S*|BNP|натрийуретическ\w*)[^.;\n]{0,40}?"
    r"(?:\bкрыс|\bмыш[еиь]|\bсвин|\bсобак|\bкош[еа]к|кошач|\bживотн)"
    r"|(?:\brat\b|\bmouse\b|\bmurine\b|porcine|bovine|canine|feline)\s+(?:IL|interleukin|NT-?proBNP|BNP)"
    r"|ветеринар|veterinar|для животных", re.I)
NOT_CONCENTRATION = re.compile(r"цитометр|cytometr|ELISpot|ЭлиСпот|иммуногистохим|immunohisto|"
                               r"ПЦР|PCR|экспресси\w+ ген|генотип|полиморфизм|амплификац", re.I)
QUANT = re.compile(r"(?<!полу)количествен|концентрац|quantitat|содержани|уровн", re.I)
QUALITATIVE = re.compile(r"качественн|полуколичествен|qualitativ|semi-?quantitat|выявлени|обнаружени", re.I)
# Cytometric Bead Array и т.п. — измеряют концентрацию цитокинов в сыворотке/плазме
CBA = re.compile(r"\bCBA\b|\bСВА\b|Capture Beads|флуоресцентных частиц|Cytometric Bead|мультиплекс|multiplex", re.I)
CALIBRATOR = re.compile(r"калибратор|calibrator|калибровочн", re.I)
CONTROL = re.compile(r"контрольн\w+ (?:материал|сыворот|образ|раствор)|контрол[ья] качества|control", re.I)
REAGENT = re.compile(r"реагент|reagent|тест|test|assay|кассет|картридж|набор|kit", re.I)


def il_numbers(text):
    nums = set()
    for m in IL_LIST.finditer(text):
        nums.update(int(n) for n in re.findall(r"\d+", m.group(1)) if len(n) <= 2)
    return nums


def markers_in(text):
    nums = il_numbers(IL_RECEPTOR.sub(" ", text))
    out = [f"IL-{n}" for n in (6, 8, 10) if n in nums]
    if NTPROBNP.search(text):
        out.append("NT-proBNP")
    return out


def relevant_lines(name, models, descr):
    """Строки состава/моделей, в которых упоминаются нужные маркеры."""
    lines = list(models) + [l.strip() for l in re.split(r"\n|<br\s*/?>", descr or "") if l.strip()]
    rel, seen = [], set()
    for l in lines:
        short = l.split(":", 1)[-1].strip() if len(l) > 250 and ":" in l else l
        if markers_in(l) and short not in seen:
            seen.add(short)
            rel.append(short)
    return rel


def classify(pid):
    card, models = details[pid]["card"], details[pid]["models"]
    name = (card.get("name") or "").strip()
    descr = card.get("modelsDescription") or ""
    func = card.get("function") or ""
    full = "\n".join([name, descr, func, *models])

    markers = markers_in(full)
    rel = relevant_lines(name, models, descr)
    # контекст для оценки образца/вида: название + релевантные строки состава
    ctx = "\n".join([name, func, *rel]) if rel else full

    if BLOOD.search(ctx):
        sample = "кровь/сыворотка/плазма"
    elif OTHER_MATRIX.search(ctx):
        sample = "другой материал (не кровь)"
    else:
        sample = "не указано"

    if CALIBRATOR.search(name) and not re.search(r"реагент", name, re.I):
        kind = "калибратор"
    elif CONTROL.search(name) and not re.search(r"реагент", name, re.I):
        kind = "контрольный материал"
    elif REAGENT.search(ctx):
        kind = "набор реагентов / тест"
    else:
        kind = "прочее"

    flags = []
    if not markers:
        flags.append("маркеры IL-6/8/10, NT-proBNP не найдены (только BNP/другие IL)")
    if ANIMAL.search(ctx):
        flags.append("не для человека (животные)")
    if CBA.search(ctx):
        flags.append("мультиплекс на частицах (CBA, проточная цитометрия) — измеряет концентрацию, "
                     "проверить тип образца")
    elif NOT_CONCENTRATION.search(ctx):
        flags.append("не иммуноанализ концентрации (цитометрия/ПЦР/ИГХ)")
    if QUALITATIVE.search(ctx) and not QUANT.search(ctx):
        flags.append("качественный/полуколичественный тест")
    if sample == "другой материал (не кровь)":
        flags.append("образец не кровь")
    if kind == "прочее":
        flags.append("тип изделия неясен")

    if not markers or any(f.startswith(("не для", "не иммуно", "образец не")) for f in flags):
        suitable = "нет"
    elif sample == "не указано" or kind == "прочее" or len(flags) > 0:
        suitable = "проверить по инструкции"
    else:
        suitable = "да"

    producer = card.get("producer") or {}
    declarant = card.get("declarant") or {}
    return {
        "Подходит": suitable,
        "Маркеры": ", ".join(markers),
        "IL-6": "IL-6" in markers,
        "IL-8": "IL-8" in markers,
        "IL-10": "IL-10" in markers,
        "NT-proBNP": "NT-proBNP" in markers,
        "Только BNP (не NT-proBNP)": (not NTPROBNP.search(full)) and bool(BNP_ONLY.search(full)),
        "Образец": sample,
        "Тип изделия": kind,
        "Количественный": bool(QUANT.search(ctx)),
        "Примечания": "; ".join(flags),
        "Номер РУ": card.get("noRu"),
        "Номер ЕРУЛ": card.get("noErul"),
        "Дата РУ": card.get("dateRu"),
        "Срок действия РУ": card.get("endDateRu") or "бессрочно",
        "Статус": (card.get("status") or {}).get("name"),
        "Наименование": name,
        "Релевантные позиции состава": "\n".join(rel),
        "Назначение": func,
        "Производитель": producer.get("name"),
        "Адрес производителя": producer.get("actualAddress") or producer.get("legalAddress"),
        "Заявитель": declarant.get("name"),
        "Места производства": "; ".join(
            f"{s.get('engName') or ''} {s.get('address') or ''}".strip()
            for s in card.get("productionSites") or []),
        "Число моделей/позиций": len(models),
        "Состав (модели)": "\n".join(models) if models else descr,
        "ID записи": pid,
        "Ссылка": f"{BASE}/widget/med-product/{pid}",
        "Найдено по запросам": ", ".join(details[pid]["terms"]),
    }


df_all = pd.DataFrame([classify(pid) for pid in details])
order = {"да": 0, "проверить по инструкции": 1, "нет": 2}
df_all = (df_all.sort_values(["Подходит", "Маркеры", "Дата РУ"],
                             key=lambda s: s.map(order) if s.name == "Подходит" else s,
                             ascending=[True, True, False])
                .reset_index(drop=True))

# Итоговая таблица: подходящие и требующие проверки
df = df_all[df_all["Подходит"] != "нет"].reset_index(drop=True)

print(df_all["Подходит"].value_counts().to_string(), "\n")
for m in ["IL-6", "IL-8", "IL-10", "NT-proBNP"]:
    print(f"{m:10} подходит/проверить: {int(df[m].sum())}")


# %% 6. Сохранение в Excel
def to_excel_safe(frame):
    # Excel ограничивает ячейку 32767 символами
    return frame.apply(lambda c: c.map(lambda v: v[:32000] if isinstance(v, str) else v))


with pd.ExcelWriter(OUTPUT_XLSX, engine="openpyxl") as xw:
    to_excel_safe(df).to_excel(xw, sheet_name="Подходящие", index=False)
    to_excel_safe(df_all).to_excel(xw, sheet_name="Все найденные", index=False)
    for m in ["IL-6", "IL-8", "IL-10", "NT-proBNP"]:
        to_excel_safe(df[df[m]]).to_excel(xw, sheet_name=m, index=False)
print("Сохранено:", os.path.abspath(OUTPUT_XLSX))

# В Colab — скачать файл:
# from google.colab import files; files.download(OUTPUT_XLSX)

df.head(20)
