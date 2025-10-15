import time
import logging
import psutil


def get_memory_usage():
    process = psutil.Process()
    return process.memory_info().rss


def format_bytes(bytes_value):
    if bytes_value < 1024 ** 2:
        return f"{bytes_value / 1024:.1f} KB"
    elif bytes_value < 1024**3:
        return f"{bytes_value / (1024**2):.1f} MB"
    else:
        return f"{bytes_value / (1024**3):.3f} GB"


def format_time(seconds):
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)

    if hours > 0:
        return f"{int(hours)}h {int(minutes)}m {int(seconds)}s"
    elif minutes > 0:
        return f"{int(minutes)}m {int(seconds)}s"
    else:
        return f"{int(seconds)}s"


def log_runtime(start_time, start_memory=None, task_name=""):
    end_time = time.time()
    end_memory = get_memory_usage() if start_memory is not None else None

    elapsed = end_time - start_time

    if task_name:
        logging.info(f"{task_name} - Runtime: {format_time(elapsed)}")
    else:
        logging.info(f"Runtime: {format_time(elapsed)}")

    if start_memory is not None and end_memory is not None:
        memory_diff = end_memory - start_memory
        memory_diff_str = format_bytes(memory_diff) if memory_diff > 0 else f"-{format_bytes(abs(memory_diff))}"
        if task_name:
            logging.info(f"{task_name} - Memory change: {memory_diff_str} (Final: {format_bytes(end_memory)})")
        else:
            logging.info(f"Memory change: {memory_diff_str} (Final: {format_bytes(end_memory)})")
