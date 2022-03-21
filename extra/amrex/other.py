
def progressBar(iterable, prefix = 'Progress:', suffix = 'Complete', decimals = 1, length = 50, fill = '█', printEnd = "\r", enumeration=False):
    """
    Call in a loop to create terminal progress bar
    This fancy function is taken from https://stackoverflow.com/a/34325723.

    Parameters:
    ------------------------------
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\\r", "\\r\\n") (Str)

    Usage example
    ----------------------------
    \# create an iterable list of items \n
    items = list(range(0, 57))

    for item in progressBar(items, prefix = 'Progress:', suffix = 'Complete', length = 50):
    \n\# Do stuff...
    
    """
    total = len(iterable)
    # Progress Bar Printing Function
    def printProgressBar (iteration):
        percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
        filledLength = int(length * iteration // total)
        bar = fill * filledLength + '-' * (length - filledLength)
        print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Initial Call
    printProgressBar(0)
    # Update Progress Bar
    for i, item in enumerate(iterable):
        if enumeration:
          yield i, item
        else:
          yield item
        printProgressBar(i + 1)
    # Print New Line on Complete
    print("\n")

def set_size(width, fraction=1, subplots=(1, 1)):
    """Set figure dimensions to avoid scaling in LaTeX.

    Parameters
    ----------
    width: float or string
            Document width in points, or string of predined document type
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy
    subplots: array-like, optional
            The number of rows and columns of subplots.
    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    if width == 'thesis':
        width_pt = 426.79135
    elif width == 'beamer':
        width_pt = 307.28987
    else:
        width_pt = width

    # Width of figure (in pts)
    fig_width_pt = width_pt * fraction
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    # https://disq.us/p/2940ij3
    golden_ratio = (5**.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio * (subplots[0] / subplots[1])

    return (fig_width_in, fig_height_in)

def import_file_as_module(module_path, module_name):
  """
    Import a regular file as python module, following the python syntax: \n
    import module_path as module_name

    Parameters:
    ------------------------------
        module_path   - Required  : path to file with filename (str)
        module_name   - Required  : name of the module (str)

    Usage example
    ----------------------------
    \# import the inputfile as module inputfile \n
    import_file_as_module(inputFilePath+'SEC_Plenum_Arrhenius.py', 'inputfile') \n
    \# now we can import some objects from this new module \n
    from inputfile import T_ref, ControlOptions \n
    print(ControlOptions)
  """
  import sys, os
  if not os.path.isfile(module_path):
   raise FileNotFoundError('given file: {} does not exist!'.format(module_path))

  import importlib.machinery
  import importlib.util
  loader = importlib.machinery.SourceFileLoader(module_name, module_path)
  spec = importlib.util.spec_from_loader(loader.name, loader)
  mod = importlib.util.module_from_spec(spec)
  loader.exec_module(mod)
  sys.modules[module_name] = mod
  return mod