import os

env = Environment(CXX='clang++', 
	CCFLAGS=['-std=c++11', '-stdlib=libc++'], 
	LIBPATH=['/usr/local/lib'], 
	CPPPATH=['/usr/local/include'])

# useC11 = False

# if useC11:
#     env = Environment(ENV={'PATH' : os.environ['PATH']})
#     env.Replace(CXX='clang++')
#     env.Append(CCFLAGS = ['-std=c++11'])

env.Append(CCFLAGS=['-O3', '-g3', '-ffast-math'], CPPPATH='#turbooctospice', LIBPATH='#turbooctospice')

# # Check for required libraries unless we're cleaning up
if not env.GetOption('clean'):
    conf = Configure(env)
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