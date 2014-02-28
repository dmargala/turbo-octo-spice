env = Environment(CCFLAGS = '-O2', CPPPATH='#turbooctospice', LIBPATH='#turbooctospice')

# Check for required libraries unless we're cleaning up
if not env.GetOption('clean'):
	conf = Configure(env)
	if not conf.CheckLib('likely'):
	    print 'Did not find likely!'
	    Exit(1)
	if not conf.CheckLib('cosmo'):
		print 'Did not find cosmo!'
		Exit(1)
	if not conf.CheckLib('bosslya'):
		print 'Did not find bosslya!'
		Exit(1)
	if not conf.CheckLib('boost_program_options'):
		print 'Did not find program options'
		Exit(1)
	env = conf.Finish()

# Use alternate build directory, do not copy files
SConscript('src/SConscript', variant_dir='build', duplicate=0, exports = ['env'])