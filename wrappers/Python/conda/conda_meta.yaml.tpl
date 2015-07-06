package:
  name: coolprop
  version: {{ version }}

source:
  fn: {{ fil }}
  url: {{ url }}
  md5: {{ md5 }}

  # If this is a new build for the same version, increment the build
  # number. If you do not include this key, it defaults to 0.
  # number: 1

requirements:
  build:
    - python
    - setuptools{% for pkg in pkgs %}
    - {{ pkg -}}
{% endfor %}

  run:
    - python{% for pkg in pkgs %}
    - {{ pkg -}}
{% endfor %}

test:
  # Python imports
  imports:
    - CoolProp
    #- CoolProp.GUI
    #- CoolProp.Plots
    - CoolProp.tests

  # commands:
    # You can put test commands to be run here.  Use this to test that the
    # entry points work.


  # You can also put a file called run_test.py in the recipe that will be run
  # at test time.

  # requires:
    # Put any additional test requirements here.  For example
    # - nose

about:
  home: {{ home }}
  license: {{ license }}
  summary: {{ sumamry }}

# See
# http://docs.continuum.io/conda/build.html for
# more information about meta.yaml