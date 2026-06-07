from __future__ import print_function


def run():
    """Run the CoolProp test suite (used by ``CoolProp.test()``).

    Migrated from nose (dead on Python 3.12+) to pytest; requires pytest to be
    installed (a test-only dependency).
    """
    import os

    try:
        import pytest
    except ImportError:
        raise ImportError("Running the CoolProp tests requires pytest (pip install pytest)")

    print('about to run the tests, please be patient')
    this_path, _file = os.path.split(os.path.abspath(__file__))
    return pytest.main([this_path])


if __name__ == '__main__':
    run()
