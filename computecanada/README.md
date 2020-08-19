# seqtools on Compute Canada servers

*All commands use Beluga server*

## Install

<a name="connect"></a>

### Connect to the server

Use SSH command on Mac, Linux and Windows 10 (PowerShell). On other versions of Windows, use [Putty](https://www.putty.org)

```
ssh beluga.computecanada.ca
```

### Run the configuration script

```
curl https://raw.githubusercontent.com/francoisrobertlab/seqtools/master/bash/configure_seqtools.sh >> configure_seqtools.sh
chmod 744 configure_seqtools.sh
./configure_seqtools.sh
```

### Reconnect to the server

*If using Putty, close the window and reconnect*

On Mac, Linux and Windows 10

```
exit
ssh beluga.computecanada.ca
```

### Run installation script

```
module load seqtools
install.sh
```

### Try seqtools

```
seqtools --help
```
