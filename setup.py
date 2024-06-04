from setuptools import setup, find_packages

setup(
    name='pankegg',
    version='0.1.0',
    packages=find_packages(include=['lib', 'lib.*']),
    include_package_data=True,
    py_modules=['pankegg_make_db', 'pankegg_app'],
    install_requires=[
        'flask',
        'pandas',
        'numpy',
        'scikit-learn',
        'scipy'
    ],
    entry_points={
        'console_scripts': [
            'pankegg_app=pankegg_app:start_server',
            'pankegg_make_db=pankegg_make_db:main',
        ],
    },
    package_data={
        'lib': ['../data/*.db', '../templates/*.html', '../data/*.csv', '../data/*.tsv', '../data/*.txt'],
    },
    author='Arnaud Vanbelle & Renaud Van Damme',
    author_email='arnaudvanbelle@live.be',
    description='Application Flask pour PanKEGG',
)
