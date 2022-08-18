#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#

use FileHandle;
use Data::Dumper;
use FindBin qw ($Bin);
use File::Find;
use File::Basename;
use lib $Bin;
use Fxtran;

use strict;

my @sys = qw (OMP_LIB CUDAFOR IEEE_ARITHMETIC);

sub slurp
{
  my $f = shift;
  do { local $/ = undef; my $fh = 'FileHandle'->new ("<$f"); <$fh> };
}

my %deps;

sub deps
{
  my ($f, $m2f) = @_;;

  my %sys = map { ($_, 1) } @sys;

  unless ($deps{$f})
    {
      my $code = &slurp ($f);
      my @inc = ($code =~ m/\n\s*(?:INCLUDE|#include)\s+"(.*?)"/goms);
      my @use = map { uc ($_) } ($code =~ m/\n\s*use\s+(\w+)/gomis);
      @use = grep { ! $sys{$_} } @use;
      my @mod;
      for my $use (@use)
        {
          if (my $f90 = $m2f->{$use})
            {
              ($f90 = &basename ($f90)) =~ s/\.F90$//o;
              push @mod, $f90 . '.o';
            }
          else
            {
              push @mod, lc ($use) . '.o';
            }
        }
      $deps{$f} = [@mod, @inc];
   }

  return @{ $deps{$f} };
}

my %kindname;

sub kindname
{
  my $f = shift;

  unless ($kindname{$f})
    {
      my @code = split (m/\n/o, &slurp ($f));
     
      my ($kind, $name);
      for (@code)
        {
          ($kind, $name) = (m/^\s*(PROGRAM|SUBROUTINE|MODULE|FUNCTION)\s+(\w+)/io);
          last if ($kind);
        }
      die $f unless ($kind);
     
      @{ $kindname{$f} } = (uc ($kind), uc ($name))
    }

  return  @{ $kindname{$f} };
}

my $fh = 'FileHandle'->new (">Makefile");

$fh->print (<< "EOF");
include Makefile.inc

EOF

my (@obj, @prg);

my @F90;

for my $dir (@ARGV)
  {
    &find ({wanted => sub { my $f = $File::Find::name; return unless ($f =~ m/\.F90$/o); push @F90, $f }, no_chdir => 1}, $dir);
  }

my %mod2F90;

for my $F90 (@F90)
  {
    my ($kind, $name) = &kindname ($F90);
    if ($kind eq 'MODULE')
      {
        $mod2F90{$name} = $F90;
      }
    elsif ($kind eq 'PROGRAM')
      {
        push @prg, $F90;
      }
  }

for my $F90 (@F90)
  {
    my ($kind, $name) = &kindname ($F90);
    if ($kind eq 'SUBROUTINE')
      {
        my $intfb = &Fxtran::intfb ($F90);
        (my $ok = $intfb) =~ s/\.h$/.ok/o;
        my @deps = &deps ($intfb, \%mod2F90);
        $fh->print (<< "EOF");
$ok: $intfb @deps
	touch $ok

EOF
      }
  }

for my $h (grep { ! m/\.intfb\.h$/o } <*.h>)
  {
    (my $ok = $h) =~ s/\.h$/.ok/o;
    $fh->print (<< "EOF");
$ok: $h
	touch $ok
EOF
  }

my @subobj;

for my $F90 (@F90)
  {
    my @deps = &deps ($F90, \%mod2F90);
    for (@deps)
      {
        s/\.h$/.ok/o;
      }
    (my $obj = $F90) =~ s/\.F90$/.o/o;
    $obj = &basename ($obj);
    $fh->print (<< "EOF");
$obj: $F90 @deps
	\$(FC) -c $F90

EOF
    my ($kind) = &kindname ($F90);
    
    push @obj, $obj if ($kind ne 'PROGRAM');

    push @subobj, $obj if ($kind ne 'MODULE');
  }

for my $prg (@prg)
  {
    (my $exe = $prg) =~ s/\.F90$/.x/o;
    $exe = &basename ($exe);
    (my $obj = $prg) =~ s/\.F90$/.o/o;
    $obj = &basename ($obj);
    $fh->print (<< "EOF");
$exe: $obj @obj
	\$(FC) -o $exe $obj @obj \$(LIBS)

EOF
  }

$fh->print (<< "EOF");

subclean:
	\\rm -f @subobj

clean:
	\\rm -f *.o *.xml *.a *.x *.mod *.optrpt

tidy:
	\\rm -f *.xml *.optrpt
EOF

