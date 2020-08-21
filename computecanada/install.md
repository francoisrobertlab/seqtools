# Install seqtools on Compute Canada servers

:memo: *The examples use Beluga server*

## Connect to the server

Use SSH command inside a terminal on [Mac](https://support.apple.com/en-ca/guide/terminal/apd5265185d-f365-44cb-8b09-71a064a42125/mac), Linux or [Windows 10 (PowerShell)](https://www.howtogeek.com/662611/9-ways-to-open-powershell-in-windows-10/)

On older versions of Windows, use [Putty](https://www.putty.org)

```
ssh beluga.computecanada.ca
```

## Run the configuration script

```
curl https://raw.githubusercontent.com/francoisrobertlab/seqtools/master/bash/configure_seqtools.sh >> configure_seqtools.sh
chmod 744 configure_seqtools.sh
./configure_seqtools.sh $email@ircm.qc.ca
```

Replace `$email@ircm.qc.ca` with your email address

## Reconnect to the server

*If using Putty, close the window and reconnect*

On Mac, Linux and Windows 10

```
exit
ssh beluga.computecanada.ca
```

## Run installation script

```
module load seqtools
install.sh
```

## Try seqtools

```
seqtools --help
```
