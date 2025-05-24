import argparse
import math
import matplotlib.pyplot as plt
import numpy as np
import logging
from pathlib import Path
from typing import List, Dict  # Добавлен импорт Dict

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    filename='zhr_analysis.log',
    encoding='utf-8'
)
logger = logging.getLogger(__name__)

def calculate_r(magnitudes: List[float]) -> float:
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
    """Чтение данных из файла"""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            return [float(line.strip()) for line in f if line.strip()]
    except Exception as e:
        logger.error(f"Ошибка чтения файла: {str(e)}")
        raise

def calculate_statistics(zhr_values: List[float], expected_zhr: float = None) -> dict:
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
        # Исправленная валидация
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
            # Рассчет популяционного индекса
            r = calculate_r(magnitudes)
            
            # Рассчет ZHR (формула-заглушка, нужно реализовать вашу формулу)
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