from __future__ import print_function


def run():
    import nose, os

    print('about to run the tests, please be patient')
    this_path, file = os.path.split(os.path.abspath(__file__))
    nose.run(argv=['--where', this_path])


if __name__ == '__main__':
    run()
