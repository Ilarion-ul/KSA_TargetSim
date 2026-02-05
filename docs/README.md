# KSA_TargetSim

## Призначення
KSA_TargetSim — мінімальний Geant4-симулятор для моделювання ланцюжка взаємодій:
**e⁻ → γ → (γ,n) → n** у двошаровій мішені.
Проєкт призначений для швидких batch-прогонів, базової візуалізації та порівняння конфігурацій мішені.

## Конфігурації мішені
Підтримуються два пресети матеріалів:
- `W-Ta`
- `U-Al`

Геометрія: **двошаровий циліндр** (substrate + coating), співвісний із віссю `z`.

## Структура директорій (скорочено)
- `app/` — основний застосунок `ksasim` (джерела, include, макроси, JSON-конфіги)
- `analysis/root/` — ROOT-макроси для післяобробки
- `cmake/` — CMake підказки
- `docs/` — документація
- `scripts/` — допоміжні shell-скрипти (env)
- `tests/` — smoke-тести (CTest)
- `results/` — вихідні дані (`logs`, `root`, `vis`)

## Збірка (Ninja)
### Залежності
- CMake ≥ 3.20
- C++17 компілятор
- Geant4 (обов’язково)
- ROOT (необов’язково, якщо потрібен ROOT-вивід)

### Підготовка середовища
```bash
source scripts/env.sh
```

### Конфігурація, збірка, тести
```bash
cmake -S . -B build -G Ninja -DKSA_USE_ROOT=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

## Запуски
```bash
# Довідка
build/bin/ksasim --help

# Batch
build/bin/ksasim -m app/macros/run_batch.mac -c app/config/default_WTa.json

# Візуалізація
build/bin/ksasim --vis -c app/config/quick_vis.json
```

## JSON-конфіг
Основні секції:
- `beam`: енергія, розміри пучка, дивергенція, стартова позиція/напрям, дефекти
- `target`: тип мішені, геометричні параметри, температура
- `run`: кількість подій, потоки, вихідні каталоги/файли, візуалізація
- `physics`: physics list, cut, майбутні перемикачі photonuclear
- `geometry`: параметри world та спрощеної геометрії
- `defects` (резерв / TODO): деталізація виробничих дефектів

Приклади полів див. у:
- `app/config/default_WTa.json`
- `app/config/default_UAl.json`
- `app/config/quick_vis.json`

## Вихідні дані
- `results/logs` — текстові підсумки (JSON)
- `results/root` — ROOT-файли (якщо ROOT доступний)
- `results/vis` — збережені зображення/графіки

Нормування: метрики слід інтерпретувати **per primary e⁻**.

## Валідація та регресія
- Smoke-тести: `tests/` (CTest + `test_smoke.cmake`)
- TODO: додати `docs/validation.md` для regression notes і контрольних чисел між релізами.

## Нотатка щодо OpenGL
Під SSH або при fallback на `llvmpipe` OpenGL-візуалізація може бути повільною.
Для прискорення вимикайте траєкторії у візуальному макросі.
