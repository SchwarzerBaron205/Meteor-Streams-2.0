# Meteor-Streams-2.0

Этот репозиторий содержит проект для расчета метеорных потоков версии 2.0.

ЧТО НОВОГО:

1. Написано на языке Python;

2. По сравнению с предыдущей версией была улучшена статистическая часть, добавлен функционал построения графика динамики. Улучшен вывод;

ZHR (Zenithal Hourly Rate - Зенитное Часовое Число) - это количество метеоров в час, которое мы могли бы увидеть, если бы радиант метеорного потока находился в зените и небо было ясным.

Программа получает файл со звездными величинами метеоров и выдаёт файл в формате:

=== ОТЧЁТ О НАБЛЮДЕНИЯХ МЕТЕОРНОГО ПОТОКА ===

ОБЩИЕ ПАРАМЕТРЫ:
1. Название файла, откуда считываются данные;
2. Количество сессий за ночь (одна сессия = 1 час);
3. Количество увиденных метеоров.

СТАТИСТИЧЕСКИЕ ПОКАЗАТЕЛИ:
1. Среднее значение ZHR за ночь;
2. Стандартное отклонение среднего ZHR от введённой заранее нормы;
3. Минимальное ZHR;
4. Максимальное ZHR;
5. Дисперсию наблюдений.

СРАВНЕНИЕ С ТЕОРИЕЙ:
1. Ожидаемое ZHR (задаётся вручную);
2. Систематическое отклонение;
3. Относительная погрешность.

ПОЧАСОВЫЕ РЕЗУЛЬТАТЫ:

Час N:
• Количество метеоров;
• Высота радианта;
• Коэф. прозрачности;
• Предельная величина;
• Популяционный индекс (r);
► ZHR;
► Отклонение от среднего.

Пример использования:

Также помимо текстового файла с результатами прилагается график динамики ZHR в формате .png.

Пример запуска:

python ZHR_RESEARCH.py astro.txt --hours 3 --counts 15 2 15 --hr 40 45 50 --k 0.9 0.85 0.8 --lm 6.5 6.3 6.0 --expected-zhr 68.5 --output report.txt --plot zhr_dynamics.png

ВАЖНО:

Если количество аргументов counts, hr, k, lm не совпадает с hours (в примере запуска: hours = 3, количество аргументов везде также 3), то код выдаст ошибку!

Сам код ниже:

"""
Модуль для анализа данных метеорных потоков и расчета часового числа метеоров (ZHR).

Основные функции:
- Расчет популяционного индекса r
- Загрузка данных наблюдений из файла
- Статистическая обработка результатов
- Генерация текстового отчета
- Визуализация динамики ZHR

Пример использования:
    python zhr_analysis.py data.txt --hours 3 --counts 15 2 15 
    --hr 40 45 50 --k 0.9 0.85 0.8 --lm 6.5 6.3 6.0 
    --expected-zhr 68.5 --output report.txt --plot graph.png
"""

import argparse
import math
import matplotlib.pyplot as plt
import numpy as np
import logging
from pathlib import Path
from typing import List, Dict

# [Остальные импорты и код без изменений...]

def calculate_r(magnitudes: List[float]) -> float:
    """Расчитывает популяционный индекс r со сглаживанием Лапласа.
    
    Args:
        magnitudes (List[float]): Список звездных величин метеоров 
                                  в диапазоне от -8 до +10
        
    Returns:
        float: Значение популяционного индекса r (обычно 2.0-3.5)
        
    Raises:
        ValueError: Если менее 3 метеоров или недостаточно данных
                    для расчета отношений гистограмм
    """
    
     """Расчет популяционного индекса r со сглаживанием"""
     
    if len(magnitudes) < 3:
        raise ValueError("Слишком мало данных для расчета")
    
    int_mags = np.round(magnitudes).astype(int)
    min_mag = int(np.min(int_mags))
    max_mag = int(np.max(int_mags))
    
    bins = np.arange(min_mag - 0.5, max_mag + 1.5)
    hist, _ = np.histogram(int_mags, bins=bins)
    hist = hist + 1  # Сглаживание Лапласа
    
    ratios = []
    for m in range(min_mag, max_mag):
        n_m = hist[m - min_mag]
        n_m1 = hist[m - min_mag + 1]
        ratios.append(n_m1 / n_m)
    
    if not ratios:
        raise ValueError("Недостаточно данных для расчета r")
    
    return np.mean(ratios)


def load_data(file_path: Path) -> List[float]:
    """Загружает данные наблюдений из текстового файла.
    
    Args:
        file_path (Path): Путь к файлу с данными в формате:
                          одна звездная величина на строку
                          
    Returns:
        List[float]: Список звездных величин метеоров
        
    Raises:
        IOError: При проблемах с чтением файла
        ValueError: При некорректных числовых значениях
    """
    """Чтение данных из файла"""
    try:
    
        with open(file_path, 'r', encoding='utf-8') as f:
        
            return [float(line.strip()) for line in f if line.strip()]
            
    except Exception as e:
    
        logger.error(f"Ошибка чтения файла: {str(e)}")
        
        raise

def calculate_statistics(zhr_values: List[float], expected_zhr: float = None) -> dict:
    """Вычисляет статистические показатели для значений ZHR.
    
    Args:
        zhr_values (List[float]): Список значений ZHR по часам
        expected_zhr (float, optional): Теоретически ожидаемое значение ZHR
        
    Returns:
        dict: Словарь статистических показателей:
            - mean: среднее арифметическое
            - std_dev: стандартное отклонение
            - min: минимальное значение
            - max: максимальное значение
            - variance: дисперсия
            - deviations: абсолютные отклонения от среднего
            - bias: абсолютное отклонение от ожидаемого ZHR (если указан)
            - relative_error: относительная погрешность в % (если указан expected_zhr)
    """
    
    """Расчёт статистических показателей"""
    
    stats = {
        'mean': np.nanmean(zhr_values),
        'std_dev': np.nanstd(zhr_values),
        'min': np.nanmin(zhr_values),
        'max': np.nanmax(zhr_values),
        'variance': np.nanvar(zhr_values),
        'deviations': [abs(zhr - np.nanmean(zhr_values)) for zhr in zhr_values]
    }
    
    if expected_zhr is not None:
    
        stats['bias'] = abs(np.nanmean(zhr_values) - expected_zhr)
        stats['relative_error'] = (stats['bias'] / expected_zhr) * 100
    
    return stats

def save_report(results: List[Dict], params: Dict, stats: Dict, output_file: Path):
    """Генерирует текстовый отчет в формате .txt.
    
    Args:
        results (List[Dict]): Почасовые результаты расчета с ключами:
                              - r: популяционный индекс
                              - zhr: часовое число метеоров
        params (Dict): Параметры наблюдений:
                      - input_file: путь к исходным данным
                      - hours: количество часов наблюдений
                      - counts: количество метеоров по часам
                      - hr: высота радианта (градусы)
                      - k: коэффициент прозрачности (0.0-1.0)
                      - lm: предельная звездная величина
                      - expected_zhr: ожидаемое значение ZHR (опционально)
        stats (Dict): Результаты calculate_statistics
        output_file (Path): Путь для сохранения отчета
        
    Raises:
        IOError: При ошибках записи в файл
    """
    
    """Сохранение отчёта с расширенной статистикой"""
    
    try:
    
        with open(output_file, 'w', encoding='utf-8') as f:
        
            f.write("=== ОТЧЁТ О НАБЛЮДЕНИЯХ МЕТЕОРНОГО ПОТОКА ===\n\n")
            
            f.write("ОБЩИЕ ПАРАМЕТРЫ:\n")
            f.write(f"Файл данных: {params['input_file']}\n")
            f.write(f"Часов наблюдения: {params['hours']}\n")
            f.write(f"Всего метеоров: {sum(params['counts'])}\n\n")
            
            f.write("СТАТИСТИЧЕСКИЕ ПОКАЗАТЕЛИ:\n")
            f.write(f"Средний ZHR: {stats['mean']:.2f}\n")
            f.write(f"Стандартное отклонение: {stats['std_dev']:.2f}\n")
            f.write(f"Минимальный ZHR: {stats['min']:.2f}\n")
            f.write(f"Максимальный ZHR: {stats['max']:.2f}\n")
            f.write(f"Дисперсия: {stats['variance']:.2f}\n")
            
            if 'bias' in stats:
            
                f.write(f"\nСРАВНЕНИЕ С ТЕОРИЕЙ:\n")
                f.write(f"Ожидаемый ZHR: {params['expected_zhr']:.2f}\n")
                f.write(f"Систематическое отклонение: {stats['bias']:.2f}\n")
                f.write(f"Относительная погрешность: {stats['relative_error']:.2f}%\n")
            
            f.write("\nПОЧАСОВЫЕ РЕЗУЛЬТАТЫ:\n")
            
            for idx, res in enumerate(results, 1):
            
                f.write(f"\nЧас {idx}:\n")
                f.write(f"• Метеоров: {params['counts'][idx-1]}\n")
                f.write(f"• Высота радианта: {params['hr'][idx-1]}°\n")
                f.write(f"• Коэф. прозрачности: {params['k'][idx-1]}\n")
                f.write(f"• Предельная величина: {params['lm'][idx-1]}\n")
                f.write(f"• Популяционный индекс (r): {res['r']:.2f}\n")
                f.write(f"► ZHR: {res['zhr']:.2f}\n")
                f.write(f"► Отклонение от среднего: {stats['deviations'][idx-1]:.2f}\n")
                
    except Exception as e:
    
        logger.error(f"Ошибка записи отчёта: {str(e)}")
        
        raise

def main():
    """Обрабатывает аргументы командной строки и управляет расчетами.
    
    Поддерживаемые аргументы:
        input_file (Path): Обязательный. Файл с данными наблюдений
        --hours (int): Обязательный. Количество часов наблюдений
        --counts (int[]): Обязательный. Количество метеоров по часам
        --hr (float[]): Обязательный. Высота радианта по часам
        --k (float[]): Обязательный. Коэффициенты прозрачности
        --lm (float[]): Обязательный. Предельные звездные величины
        --expected-zhr (float): Опционально. Ожидаемое значение ZHR
        --output (Path): Обязательный. Файл для отчета
        --plot (Path): Опционально. Файл для сохранения графика
        
    Raises:
        ValueError: При несоответствии количества параметров
    """
    parser = argparse.ArgumentParser(
        description='Расчёт ZHR с расширенной статистикой',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('input_file', type=Path, help='Файл с данными')
    parser.add_argument('--hours', type=int, required=True, help='Количество часов')
    parser.add_argument('--counts', type=int, nargs='+', required=True, help='Метеоров по часам')
    parser.add_argument('--hr', type=float, nargs='+', required=True, help='Высота радианта')
    parser.add_argument('--k', type=float, nargs='+', required=True, help='Коэф. прозрачности')
    parser.add_argument('--lm', type=float, nargs='+', required=True, help='Предельная величина')
    
    parser.add_argument('--expected-zhr', type=float, help='Ожидаемое значение ZHR для сравнения')
    parser.add_argument('--output', type=Path, required=True, help='Файл отчёта')
    parser.add_argument('--plot', type=Path, help='График динамики')

    args = parser.parse_args()

    try:
        
        params_lens = {
            len(args.counts), 
            len(args.hr), 
            len(args.k), 
            len(args.lm)
        }
        
        if len(params_lens) != 1 or args.hours != len(args.counts):
            raise ValueError("Несоответствие количества параметров")

        magnitudes = load_data(args.input_file)
        results = []
        
        # Расчет ZHR для каждого часа
        for hour in range(args.hours):
        
            # Расчет популяционного индекса
            r = calculate_r(magnitudes)
            
            # Расчет ZHR = (Nkr^(6.5 - lm))/(sin(hr))
            sin_hr = math.sin(math.radians(args.hr[hour]))
            zhr = (args.counts[hour] * args.k[hour] * r**(6.5 - args.lm[hour])) / (sin_hr)
            
            results.append({'r': r, 'zhr': zhr})

        zhr_values = [res['zhr'] for res in results]
        stats = calculate_statistics(zhr_values, args.expected_zhr)

        params = {
            'input_file': args.input_file,
            'hours': args.hours,
            'counts': args.counts,
            'hr': args.hr,
            'k': args.k,
            'lm': args.lm,
            'expected_zhr': args.expected_zhr
        }
        
        save_report(results, params, stats, args.output)

        if args.plot:
            plt.figure(figsize=(10, 5))
            hours = range(1, args.hours + 1)
            plt.plot(hours, zhr_values, 'o-')
            plt.title("Динамика ZHR")
            plt.xlabel("Час наблюдения")
            plt.ylabel("ZHR")
            plt.grid(True)
            plt.savefig(args.plot, bbox_inches='tight')

    except Exception as e:
        logger.error(f"Фатальная ошибка: {str(e)}")
        print(f"Ошибка: {str(e)} (см. лог-файл)")

if __name__ == "__main__":
    main()
