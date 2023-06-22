from setuptools import find_packages, setup

setup(
    name='pyclone-vi',
    version='0.1.3',
    description='Fast method for inferring clonal population structure from SNV data.',
    author='Andrew Roth',
    author_email='aroth@bccrc.ca',
    url='https://github.com/Roth-Lab/pyclone-vi',
    packages=find_packages(),
    license='GPL v3',
    entry_points={
        'console_scripts': [
            'pyclone-vi= pyclone_vi.cli:main',
        ]
    }
)
