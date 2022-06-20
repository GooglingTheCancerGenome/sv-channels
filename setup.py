from setuptools import setup

setup(
    name='svchannels',
    version='0.0.1',
    entry_points = {
        'console_scripts': ['svchannels=svchannels.main:main'],
    },
    packages=['svchannels'],
    zip_safe=False,
)
