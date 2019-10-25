from setuptools import setup, find_packages

setup(
    name="SeqTools",
    version="0.3-SNAPSHOT",
    packages=find_packages(),
    author="Christian Poitras",
    author_email="christian.poitras@ircm.qc.ca",
    description="Tools to analyze NGS data",
    keywords="bioinformatics, NGS",
    url="https://github.com/francoisrobertlab/seqtools",
    license='GNU General Public License version 3',
    classifiers=[
        'License :: OSI Approved :: GNU General Public License version 3'
    ],
    install_requires=[
        'click>=7.0',
        'pandas>=0.25.0'
    ],
    entry_points={
        'console_scripts': [
            'bam2bed = seqtools.Bam2Bed:main',
            'download = seqtools.DownloadSample:main',
            'genecov = seqtools.GenomeCoverage:main',
            'merge = seqtools.MergeBed:main',
            'plot2do = seqtools.Plot2do:main',
            'runbwa = seqtools.RunBwa:main',
            'runvap = seqtools.RunVap:main',
            'split = seqtools.SplitBed:main',
            'statistics = seqtools.StatisticsFile:main',
        ]
    }
)
