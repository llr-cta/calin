# calin #

__Calin__ depends on certain packages to provide functionality. Some of these must be installed on the host system in order to compile and run __calin__. We refer to these as **external dependencies**. Some external dependencies may be optional. Certain functionality may be missing if they are not present, but in the build will proceed without them. These are referred to as **optional dependencies**.

Other dependencies are shipped with the __calin__ distribution in order to make building and installation easier on the user. We refer to these as **packaged dependencies**. Generally these will be smaller packages, and they must have licenses compatible with GPLv2, see (http://www.gnu.org/licenses/license-list.en.html) for details.
 
## External dependencies ##

### Google protobuf library ###

**Version:** 3.0.0-beta1

### Python and Numpy ###

### GNU scientific library (GSL) ###

## Optional dependencies ##

### Geant4 ###

## Packaged dependencies ##

### Eigen ###

**Version:** 3.2.7

**License:** MPL2 -  compatible with GPLv2

**URL:** http://eigen.tuxfamily.org/

    $ bzcat eigen-3.2.6.tar.bz2 | tar xf -
    $ cd eigen-3.2.6
	$ mv Eigen ~/calin/include
	$ mv unsupported ~/calin/include/Eigen

### NLOpt ###

### F2C ###

### Minuit-76 ###

### CMinpack ###

### GTest ###
