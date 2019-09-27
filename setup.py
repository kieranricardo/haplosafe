from setuptools import setup

REQUIRES = [
    "numpy>=1.14",
    "networkx>=2.2"
]


setup(
    name='haplosafe',
    version='0.1',
    packages=['haplosafe'],
    url='',
    license='',
    author='kieranricardo',
    author_email='',
    description='',
    python_requires=">=3.5.0, <3.8.0",
    install_requires=REQUIRES
)