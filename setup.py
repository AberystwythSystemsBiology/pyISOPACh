from distutils.core import setup

setup(
    name='pyIDICk',
    version='0.0.1',
    packages=['pyidick'],
    url='https://www.github.com/pyiDICk',
    license='MIT',
    author='Keiron OShea',
    author_email = 'keo7@aber.ac.uk',
    description = 'Isotopic Distribution Calculator',
    package_data={
      'pyidick': ['data/*.json'],
   },
)
