from setuptools import setup

setup(name='tracking_forces',
      version=0.1,
      description='compute reaction forces from tracked whisker images',
      author='Elliot Hevel',
      author_email='elliothevel2013@u.northwestern.edu',
      package_dir={'tracking_forces': 'src'},
      packages=['tracking_forces'],
      package_data={'':['tracking_forces.cfg']},
      entry_points={'console_scripts': ['tracking_forces=tracking_forces.tracking_forces:main']}
     )
