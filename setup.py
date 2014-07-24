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
      scripts=['scripts/f52fasta.py'],
      entry_points = {
        'console_scripts': [
            'coverageFromBlast6 = jptools.coverage:coverage_from_blast6',
            'correctOxford = jptools.correct:correct_oxford'
            ]
        },
      )
      
      
      
      
