"""
Entrypoint module, in case you use `python -m click_mergevcfs`.

Why does this file exist, and why __main__? For more info, read:

- https://www.python.org/dev/peps/pep-0338/
- https://docs.python.org/3/using/cmdline.html#cmdoption-m
"""

from click_mergevcfs.cli import main

if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
