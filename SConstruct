# SConstruct is python -*-Python-*-
import os
import distutils.sysconfig

# initialize our build enviroment
env = Environment(CCFLAGS=['-std=c++11'], CPPPATH=['/usr/local/include'],LIBPATH=['/usr/local/lib'])

# take PATH and LD_LIBRARY_PATH from the environment so currently configured
# build tools are used
if 'PATH' in os.environ:
   env['ENV']['PATH'] = os.environ['PATH']
if 'LD_LIBRARY_PATH' in os.environ:
   env['ENV']['LD_LIBRARY_PATH'] = os.environ['LD_LIBRARY_PATH']

# take CC and CXX from the current environment if they are defined
if 'CC' in os.environ:
   env['CC'] = os.environ['CC']
   print 'Using CC =',env['CC']
if 'CXX' in os.environ:
   env['CXX'] = os.environ['CXX']
   print 'Using CXX =',env['CXX']

# configure the environment to find the packages we need
conf = Configure(env)

# copy any library paths defined in $LDFLAGS to LIBPATH
if 'LDFLAGS' in os.environ:
   for token in os.environ['LDFLAGS'].split():
      if token[:2] == '-L':
         conf.env.Append(LIBPATH=[token[2:]])

env.Append(CCFLAGS=['-O3', '-g3', '-ffast-math'], CPPPATH='#turbooctospice', LIBPATH='#turbooctospice')

# # Check for required libraries unless we're cleaning up
if not env.GetOption('clean'):
    if not conf.CheckLib('hdf5', language='cxx'):
    	print 'Did not find hdf5!'
    	Exit(1)
    if not conf.CheckLibWithHeader('hdf5_cpp', 'H5Cpp.h', 'cxx'):
    	print 'Did not find hdf5_cpp!'
    	Exit(1)
    if not conf.CheckLibWithHeader('likely', 'likely/likely.h', 'cxx'):
        print 'Did not find likely!'
        Exit(1)
    if not conf.CheckLibWithHeader('cosmo', 'cosmo/cosmo.h', 'cxx'):
        print 'Did not find cosmo!'
        Exit(1)
    # if not conf.CheckLib('bosslya'):
    #     print 'Did not find bosslya!'
    #     Exit(1)
    if not conf.CheckLib('boost_program_options', language='cxx'):
        print 'Did not find boost_program_options'
        Exit(1)
    if not conf.CheckLib('boost_system', language='cxx'):
        print 'Did not find boost_system'
        Exit(1)
    if not conf.CheckLib('boost_thread', language='cxx'):
        print 'Did not find boost_thread'
        Exit(1)
    if not conf.CheckLib('boost_coroutine', language='cxx'):
        print 'Did not find boost_coroutine'
        Exit(1)
    env = conf.Finish()

# Use alternate build directory, do not copy files
SConscript('src/SConscript', variant_dir='build', duplicate=0, exports=['env'])
