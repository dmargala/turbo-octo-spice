import os

env = Environment(CXX='clang++', CCFLAGS=['-std=c++11','-stdlib=libc++'])

useC11 = False

if useC11:
    env = Environment(ENV={'PATH' : os.environ['PATH']})
    env.Replace(CXX='clang++')
    env.Append(CCFLAGS = ['-std=c++11'])

env.Append(CCFLAGS=['-O2','-g3'], CPPPATH='#turbooctospice', LIBPATH='#turbooctospice')

# # Check for required libraries unless we're cleaning up
# if not env.GetOption('clean'):
#     conf = Configure(env)
#     if not conf.CheckLib('likely'):
#         print 'Did not find likely!'
#         Exit(1)
#     if not conf.CheckLib('cosmo'):
#         print 'Did not find cosmo!'
#         Exit(1)
#     if not conf.CheckLib('bosslya'):
#         print 'Did not find bosslya!'
#         Exit(1)
#     if not conf.CheckLib('boost_program_options'):
#         print 'Did not find program options'
#         Exit(1)
#     env = conf.Finish()

# Use alternate build directory, do not copy files
SConscript('src/SConscript', variant_dir='build', duplicate=0, exports = ['env'])