# KSA_TargetSim

## Назначение и обзор
KSA_TargetSim — Geant4‑симулятор для моделирования цепочки взаимодействий
**e⁻ → γ → (γ,n) → n** в мишени. Проект ориентирован на быстрые batch‑прогоны,
сравнение конфигураций мишени и получение воспроизводимых метрик для:
- нейтронного поля (спектр, угловое распределение, карты потоков);
- энерговыделения (по объёмам, сеткам, 3D‑карте);
- подготовительных оценок DPA и газообразования (H/He).

Нормирование всех метрик — **per primary e⁻**. Масштабирование «на секунду»
выполняется в постобработке через `N_e_per_s = I_avg / e`.

## Конфигурации мишени и геометрии
Поддерживаются три типа:
- `W-Ta` — секционная мишень с пластинами W, буфером Ti, покрытием Ta и водяными зазорами.
- `U-Al` — упрощённый цилиндр (U + Al).
- `U-Mo` — реалистичная мишень с пластинами U7Mo, Al‑покрытием, водяными зазорами,
  входным окном, корпусом и гелиевой камерой.

### Имена объёмов (ключевые для скоринга)
Используются в `SteppingAction` для агрегации:
- W‑Ta: `TargetSubstrate`, `TargetCoating`, `TargetBufferTi`, `TargetWaterGap`, `TargetAssembly`.
- U‑Mo: `PlateU_<i>`, `PlateCladAl_<i>`, `WaterGap_<i>`, `EntranceWindow`, `TargetHousing`, `HeliumChamber`.

## Физическая модель (детальный блок)
### Физические процессы и список
- Базовая физика Geant4 (стандартный список, заданный в `physics.physicsListName`).
- Поддержка **photonuclear** процессов (переключатель `enablePhotonuclear`).
- Для нейтронов используется HP‑модель (в составе `QGSP_BIC_HPT`).

### Механизм генерации нейтронов
Первичный электрон генерирует γ‑каскад в мишени; далее включаются (γ,n) реакции,
в результате чего появляются вторичные нейтроны. Отслеживаются:
- энергия нейтрона `En` (MeV),
- направление (через `cosθ` относительно оси +Z),
- координаты выхода на поверхности мишени.

### Энерговыделение
Суммируется:
- по слоям (core / clad / buffer / water / housing);
- по пластинам (per‑plate);
- в 3D‑сетке (грид) для построения карт энергоплотности.

### DPA proxy (NIEL)
Используется `step->GetNonIonizingEnergyDeposit()` как **proxy** повреждений.
Эта величина требует последующей интерпретации (NRT/arc‑dpa) с порогом смещения `E_d`.
В проекте сохраняется по пластинам как “NIEL proxy”.

### Газообразование (H/He)
Счётчик вторичных частиц:
- H‑группа: `proton`, `deuteron`, `triton`;
- He‑группа: `alpha`, `He3`.

Данные агрегируются per‑plate и используются как заготовка для оценки окрихчения.

## Выходные данные
### Текстовые
- `results/logs/run_summary.json` — сводка запуска.
- `results/logs/heatmaps.json` — массивы тепловых/нейтронных карт (для не‑ROOT режима).

### ROOT
`results/root/<run_id>.root` содержит:
- `RunMeta` — метаданные запуска (параметры пучка/геометрии + нормирование).
- `run_summary` — интегральные метрики.
- `NeutronSurf` — пересечения нейтронов с поверхностями мишени.
- `edep_3d` (TH3D) — 3D карта энерговыделения (по Z ограничена стеком пластин мишени).
- `h2_edep_xy_mid` (TH2D) — 2D срез по середине Z.
- `h2_neutron_exit_*` — карты выходов нейтронов:
  - `h2_neutron_exit_xy_downstream`, `h2_neutron_exit_xy_upstream`
  - `h2_neutron_exit_yz_side_x`, `h2_neutron_exit_xz_side_y`
  - `h2_neutron_exit_side_surface` (объединённая боковая поверхность)
- `h1_niel_plate` — NIEL‑proxy per‑plate
- `h1_gas_h_plate`, `h1_gas_he_plate` — газообразование per‑plate

### Экспорт артефактов из ROOT
Скрипт выгружает PNG‑изображения всех `TH1/TH2/TH3`, табличные данные из `NeutronSurf`
и отдельный экспорт данных источника нейтронов в отдельную папку с отметкой времени
(взято из времени создания ROOT‑файла, либо из mtime файла).

```bash
scripts/export_root_artifacts.sh results/root/run_WTa.root
```

Результат: `results/root/png/<timestamp>_<target_type>/...`


Дополнительно для источника нейтронов создаются:
- `neutron_source.csv` — табличные данные (`event_id, En_MeV, x_mm, y_mm, z_mm, cosTheta, weight, time_ns, surface_id`)
- `neutron_source_spectrum_linear.png` — спектр нейтронов по энергии в линейной шкале
- `neutron_source_spectrum_log.png` — спектр нейтронов по энергии в логарифмической шкале

## Структура директорий (кратко)
- `app/` — приложение `ksasim` (код, include, макросы, JSON‑конфиги)
- `analysis/root/` — ROOT‑макросы постобработки
- `cmake/` — CMake‑подсказки
- `docs/` — документация
- `scripts/` — вспомогательные скрипты
- `tests/` — smoke‑тесты (CTest)
- `results/` — выходные данные (`logs`, `root`, `vis`)

## Сборка (Ninja)
### Зависимости
- CMake ≥ 3.20
- C++17 компилятор
- Geant4 (обязательно)
- ROOT (опционально, для ROOT‑вывода)

### Подготовка окружения
```bash
source scripts/env.sh
```

### Конфигурация, сборка, тесты
```bash
cmake -S . -B build -G Ninja -DKSA_USE_ROOT=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

## Запуски
```bash
# Помощь
build/bin/ksasim --help

# Batch
build/bin/ksasim -m app/macros/run_batch.mac -c app/config/default_WTa.json

# Визуализация
build/bin/ksasim --vis -c app/config/quick_vis.json
```

## JSON‑конфиг
Основные секции:
- `beam`: энергия, размеры пучка, дивергенция, позиция/направление
- `target`: тип мишени, геометрия, температура
- `run`: количество событий, потоки, выходные каталоги/файлы
- `physics`: physics list, cut, photonuclear
- `geometry`: параметры мира и сборки

Примеры:
- `app/config/default_WTa.json`
- `app/config/default_UAl.json`
- `app/config/default_UMo.json`
- `app/config/quick_vis.json`

## Валидация и регрессия
- Smoke‑тесты: `tests/` (CTest + `test_smoke.cmake`)
- TODO: добавить `docs/validation.md` для контрольных чисел и regression‑логов.

## Примечание по OpenGL
Под SSH или при fallback на `llvmpipe` OpenGL‑визуализация может быть медленной.
Для ускорения отключайте траектории в визуальных макросах.
