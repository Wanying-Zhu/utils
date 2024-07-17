def progress_bar(total, progress):
    '''
    Prameters:
    - total: total number of operations to be done
    - progress: number of operations completed
    '''
    percent = 100 * (progress/total)
    bar = '=' * int(percent) + '-' * int(100 - percent)
    print(f'|{bar}| {percent:.2f}%', end='\r')