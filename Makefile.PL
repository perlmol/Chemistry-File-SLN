use 5.006;
use ExtUtils::MakeMaker;

WriteMakefile(
    'NAME'		=> 'Chemistry::File::SLN',
    'VERSION_FROM'	=> 'lib/Chemistry/File/SLN.pm', # finds $VERSION
    'PREREQ_PM'		=> {qw(
        Chemistry::Mol              0.24
    )},
);
