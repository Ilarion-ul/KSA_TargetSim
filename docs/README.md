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
- U‑Mo: `PlateU_<i>`, `PlateCladAl_<i>`, `WaterGap_<i>`, `EntranceWindow`, `TargetHousing`, `HeliumChamber`, `AlignmentPin` (optional).

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


### Матеріали установки (додатково)
- **SAV-1 (САВ-1)** введено як окремий матеріал `SAV1` у геометрії (корпус, вхідне вікно, оболонки):
  - густина: `2719 кг/м3`;
  - довідкові теплофізичні константи в моделі: `Cp = 871 Дж/(кг·К)`, `k = 202.4 Вт/(м·К)`.
- **Демінералізована вода (H2O, 25°C)** введена як `KSA_DEMIN_WATER_25C` для теплоносія/уповільнювача:
  - густина: `0.997 г/см3`;
  - тиск у матеріалі: `2.69 бар`;
  - `Cp = 4184 Дж/(кг·К)`.

> Примітка: залежність густини води від температури та параметри пари (для аварійних режимів) наразі задокументовані на рівні ТЗ і можуть бути додані окремою температурно-залежною моделлю матеріалів.

### Актуальні параметри геометрії (ТЗ)
- **W-Ta**
  - 7 пластин, `plate_xy_mm` = 65.8 (допускається 66.0)
  - `plate_thicknesses_mm` = `[2.5, 2.5, 2.5, 3.5, 3.5, 5.5, 9.5]`
  - `clad_ta_mm` у діапазоні `[0.25, 0.27]`
  - `buffer_ti_mm` у діапазоні `[0.03, 0.06]`
  - `water_gap_mm` пресетно `2.0` або `1.75`
  - `w_substrate_material` = `pure_W` або `W_Fe_Ni` (1.5% Fe, 3% Ni, решта W)
- **U-Mo**
  - 12 пластин (фінальна нейтроногенеруюча конфігурація), `plate_xy_mm` = 64.0
  - `plate_thicknesses_mm` = `[2.5,2.5,2.5,2.5,3,3,4,5,7,10,14,22.5]`
  - `clad_thickness_front_mm` = 0.95 (для пластин 1–4), `clad_thickness_rest_mm` = 0.7 (для 5–12)
  - `gap_inout_mm` = 1.0
  - `inter_plate_gaps_mm` = 11 значень по `1.75` мм
  - 6-пластинний варіант підтримується як legacy/backward-compatible
- **Загальна збірка**
  - `housing_inner_xy_mm` = 66, `housing_wall_mm` = 2
  - `total_assembly_len_mm` = 2620
  - `beamline_vacuum_len_mm` = 2210
  - `entrance_window_mm` = 2
  - `helium_chamber_len_mm` = 237
  - опціональний центруючий "палець": `geometry.enable_alignment_pin` + `alignment_pin_*`

## Выходные данные
### Текстовые
- `results/logs/run_summary.json` — сводка запуска.
- `results/logs/heatmaps.json` — массивы тепловых/нейтронных карт (для не‑ROOT режима).

### ROOT
`results/root/<run_id>.root` содержит:
- `RunMeta` — метаданные запуска (параметры пучка/геометрии + нормирование).
- `run_summary` — интегральные метрики (включно з `nNeutronModelExit` — окремий лічильник нейтронів, що перетнули межу світу/моделі).
- `NeutronSurf` — пересечения нейтронов с поверхностями мишени.
- `PhotonSurf` — пересечения фотонов с поверхностями мишени.
- `edep_3d` (TH3D) — 3D карта энерговыделения (по Z ограничена стеком пластин мишени).
- `h2_edep_xy_mid` (TH2D) — 2D срез по середине Z.
- `h2_neutron_exit_*` — карты выходов нейтронов:
  - `h2_neutron_exit_xy_downstream`, `h2_neutron_exit_xy_upstream`
  - `h2_neutron_exit_yz_side_x`, `h2_neutron_exit_xz_side_y`
  - `h2_neutron_exit_side_surface` (объединённая боковая поверхность)
- `h1_niel_plate` — NIEL‑proxy per‑plate
- `h1_gas_h_plate`, `h1_gas_he_plate` — газообразование per‑plate
- `MeshData` — воксельные данные для пост-расчёта разбухания (mesh_id/voxel_id, геометрическая привязка, edep, damage proxy, flux-группы, H/He).

### Экспорт артефактов из ROOT
Скрипт выгружает PNG‑изображения всех `TH1/TH2/TH3`, табличные данные из `NeutronSurf`
и отдельный экспорт данных источника нейтронов в отдельную папку с отметкой времени
(взято из времени создания ROOT‑файла, либо из mtime файла).

```bash
scripts/export_root_artifacts.sh results/root/run_WTa.root
```

Результат: `results/root/png/<timestamp>_<target_type>/...`


Дополнительно создаются артефакты по источникам:
- `neutron_source.csv` — табличные данные (`event_id, En_MeV, x_mm, y_mm, z_mm, cosTheta, weight, time_ns, surface_id`)
- `neutron_source_spectrum_linear.png` — спектр нейтронов в диапазоне **0–5 MeV** (линейная шкала)
- `neutron_source_spectrum_log.png` — спектр нейтронов в диапазоне **2.5e-9–5 MeV** (0.0025 eV нижняя граница) (логарифмическая шкала, для теплового хвоста)
- `photon_source_spectrum_linear.png` — спектр фотонов в диапазоне **0–100 MeV** (линейная шкала)
- `photon_source_spectrum_log.png` — спектр фотонов в диапазоне **~0–100 MeV** (логарифмическая шкала, с ненулевым минимумом по оси X)
- `photon_source_spectrum_4p5_30.png` — фокусний спектр фотонів у вікні **4.5–30 MeV**
- `particle_yields_per_electron.json` — количества `photons/electron` и `neutrons/electron` (также weighted-вариант); добавлены поля `*_from_run_summary` для сверки с интегральной статистикой run_summary, включая `neutrons_model_exit_per_electron_from_run_summary`

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

# Production-oriented presets
build/bin/ksasim -m app/macros/run_batch.mac -c app/config/production_WTa.json
build/bin/ksasim -m app/macros/run_batch.mac -c app/config/production_UMo.json

# Визуализация
build/bin/ksasim --vis -c app/config/quick_vis.json
```

## JSON‑конфиг
Основные секции:
- `beam`: энергия, размеры пучка, дивергенция, позиция/направление
- `target`: тип мишени, геометрия, температура
- `run`: количество событий, потоки, выходные каталоги/файлы (`irradiation_time_s`, `enableSwellingOutput`)
- `physics`: physics list, cut, photonuclear
- `geometry`: параметры мира и сборки

Примеры:
- `app/config/default_WTa.json`
- `app/config/default_UAl.json`
- `app/config/default_UMo.json`
- `app/config/production_WTa.json`
- `app/config/production_UMo.json`
- `app/config/quick_vis.json`

## Валидация и регрессия
- Smoke‑тесты: `tests/` (CTest + `test_smoke.cmake`)
- TODO: добавить `docs/validation.md` для контрольных чисел и regression‑логов.

## Примечание по OpenGL
Под SSH или при fallback на `llvmpipe` OpenGL‑визуализация может быть медленной.
Для ускорения отключайте траектории в визуальных макросах.


## Интерфейс данных для расчёта радиационного разбухания
Новые файлы (backward-compatible, существующие `run_summary.json` и `heatmaps.json` не меняются):
- `results/logs/run_meta.json`
- `results/logs/mesh_definition.json`
- `results/logs/mesh_data.csv`

### Что содержит `mesh_definition.json`
- `voxel_namespace` з глобально унікальним `voxel_id` діапазоном
- `meshes[]` (поточний run може містити 1..N mesh-областей; зараз заповнюється `target_plate_stack`)

### Что содержит `mesh_data.csv`
Одна строка = один воксель:
- идентификаторы: `mesh_id`, `voxel_id`, `ix`, `iy`, `iz`
- координаты/объём: `x_mm`, `y_mm`, `z_mm`, `volume_cm3`
- привязка к физике/слоям: `material_id/material_name`, `volume_id/volume_name`, `layer_id/layer_name`
- поля: `edep_MeV_per_primary`, `edep_J_per_primary`, `damage_energy_eV_per_primary`, `n_events_scored`
- спектральные колонки: `flux_n_total_per_cm2_per_primary`, `flux_n_g0..flux_n_g19`
- газ: `h_prod_per_primary`, `he_prod_per_primary`

### Нормировка и пересчёт
В `run_meta.json` сохраняются поля:
- `normalization_mode = per_primary`
- `beam_current_A`, `irradiation_time_s`, `electrons_per_second`, `electrons_total`

Используйте:
- `N_e_per_s = I / e`
- `N_e_campaign = (I / e) * t`
- `Q_s(x) = Q_pp(x) * N_e_per_s`
- `Q_campaign(x) = Q_pp(x) * N_e_campaign`

### Статус физической полноты
Текущее поле `damage_energy_eV_per_primary` помечено как proxy (`derived_from=damage_proxy`),
пока не подключён отдельный transport-consistent scorer для displacement damage / DPA.
`flux_n_g*` и `h/he` в этой версии экспортируются как interface placeholders (нулевые значения),
чтобы стабилизировать downstream schema и миграцию postprocessing-кода.

Для метаданных групп нейтронов сохраняются две сетки: логарифмическая и линейная (`neutron_energy_group_edges_MeV` и `neutron_energy_group_edges_linear_MeV`).

При `enableSwellingOutput=false` mesh-интерфейс (`run_meta.json`/`mesh_definition.json`/`mesh_data.csv` и ROOT `MeshData`) не формируется.


Примечание: `nGamma` в `run_summary` и ROOT summary считается как число фотонов с энергией > 5 MeV (порог для фотоядерной релевантности).


### Як отримати mesh-дані для розрахунку розбухання
1. У конфігу запуску встановіть `run.enableSwellingOutput = true`.
2. Запустіть розрахунок (`ksasim -m app/macros/run_batch.mac -c <config.json>`).
3. Після завершення перевірте файли в `results/.../logs/`:
   - `run_meta.json`
   - `mesh_definition.json`
   - `mesh_data.csv`
4. Якщо зібрано з ROOT (`KSA_USE_ROOT=ON`), додатково буде `MeshData` у `results/.../root/*.root`.
