from setuptools import setup

def readme():
    with open('README') as f:
        return f.read()



setup(name='jptools',
      version='0.01',
      description='For Oxford',
      long_description=readme(),
      url='https://github.com/jgurtowski/jptools',
      author='James Gurtowski',
      author_email='gurtowsk@cshl.edu',
      license='GPL',
      packages=['jptools'],
      install_requires = [
        "h5py >= 2.3.1",
        "jbio >= 0.1", 
        "pbcore >= 0.6.3",
        "pbtools.pbdagcon >= 0.2.3",
        ],
      dependency_links = [
        "git+https://github.com/jgurtowski/jbio#egg=jbio-0.1",
        "git+https://github.com/PacificBiosciences/pbcore#egg=pbcore-0.6.3",
        "git+https://github.com/PacificBiosciences/pbdagcon#egg=pbtools.pbdagcon-0.2.3"
        ],
      entry_points = {
        'console_scripts': [
            'coverageFromBlast6 = jptools.coverage:coverage_from_blast6',
            'correctOxford = jptools.correct:correct_oxford',
            'blast6Filter = jptools.blast:blast6filter_main',
            'fast5ToFasta = jptools.fast5:fast5ToFasta_main',
            ]
        },
      )
      
      
      
      
