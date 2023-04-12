from functools import wraps
import time

def timeit(func):
    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        print(f'Function {func.__name__} took {total_time:.4f}s')
        return result
    return timeit_wrapper


@timeit
def calculate_something(num):
    """
    Simple function that returns sum of all numbers up to the square of num. To demonstrate the use of timeit decorator
    """
    total = sum([x for x in range(0, num**2)])
    return total

if __name__ == '__main__':
    calculate_something(10)
    calculate_something(100)
    calculate_something(10000)