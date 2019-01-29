from version import version, required_versions
from setuptools import find_packages, setup


kwargs = dict(
    name='toil-narle',
    version=version,
    description="UCSC CGL Nanopore Toil pipeiline",
    author='UCSC Computational Genomics Lab',
    author_email='tpesout@ucsc.edu',
    url="https://github.com/",
    install_requires=[x + y for x, y in required_versions.iteritems()],
    tests_require=['pytest==2.8.3'],
    package_dir={'': 'src'},
    packages=find_packages('src'),
    entry_points={
        'console_scripts': ['toil-narle = narle.narle_pipeline:main']})


setup(**kwargs)

