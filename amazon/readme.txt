

https://wiki.x2go.org/doku.php/doc:usage:x2goclient

https://datawookie.netlify.com/blog/2017/08/remote-desktop-on-an-ubuntu-ec2-instance/

https://websiteforstudents.com/connect-to-ubuntu-16-04-17-10-18-04-desktop-via-remote-desktop-connection-rdp-with-xrdp/
https://stackoverflow.com/questions/25657596/how-to-set-up-gui-on-amazon-ec2-ubuntu-server






https://stackoverflow.com/questions/25657596/how-to-set-up-gui-on-amazon-ec2-ubuntu-server

Remote Desktop on an Ubuntu EC2 Instance
2017-08-08  AWS Andrew B. Collier

A couple of options for remote access to desktop applications on a EC2 host.

Option 1: SSH with X11 Forwarding
Connect to the EC2 host using SSH with X11 forwarding enabled.

ssh -X 13.57.185.127
In the resulting session you should find that the DISPLAY environment variable is set.

echo $DISPLAY
With this in place you can launch an application on the remote host and it will show up on your local desktop. Try starting gvim (assuming that you have it installed).

Option 2: Remote Desktop
Connect via SSH.
Install a few packages.

sudo apt update
sudo apt install -y ubuntu-desktop xrdp
Edit the RDP configuration file, /etc/xrdp/xrdp.ini, on the host. Note the entry for port, which will be important for making a connection. A minimal configuration might look like this:

[globals]
bitmap_cache=yes
bitmap_compression=yes
port=3389
crypt_level=low
channel_code=1
max_bpp=24

[xrdp1]
name=sesman-Xvnc
lib=libvnc.so
username=ask
password=ask
ip=127.0.0.1
port=ask-1
In the AWS Dashboard edit the Security Group for the EC2 instance and allow inbound TCP connections on port 3389.

Restart RDP.

sudo service xrdp restart
Choose the Window Manager for RDP connections. This involves changing the contents of a user’s .xsession file. You can copy the modified .xsession into /etc/skel/ so that it will be propagated into any newly created accounts. However, you’ll need to copy it manually into existing accounts.

Select one of the Window Manager options below (there are certainly other options too!).

XFCE
sudo apt install -y xfce4 xfce4-goodies
echo xfce4-session >~/.xsession
Unity
echo unity >~/.xsession
You’re ready to connect!

On a Linux machine, connect using vinagre. You’ll need to specify the IP address for the EC2 host and the RDP port.


A connection will be initiated and you’ll be prompted to provide your password. Leave the port unchanged.



Once you’re authenticated you should see your desktop.




This can be done. Following are the steps to setup the GUI

Create new user with password login
sudo useradd -m awsgui
sudo passwd awsgui
sudo usermod -aG admin awsgui

sudo vim /etc/ssh/sshd_config # edit line "PasswordAuthentication" to yes

sudo /etc/init.d/ssh restart


history of commands run to set up server


    1  sudo apt install dos2unix
    2  sudo bash amazon-ec2-hermes-1.sh 
    3  docker
    4  pwd
    5  ls
    6  cp amazon-ec2-hermes-1.sh amazon-ec2-hermes-2.sh 
    7  sudo bash amazon-ec2-hermes-2.sh 
    8  ls
    9  sudo /bin/rm -rf proj
   10  bash amazon-ec2-hermes-2.sh 
   11  spack
   12  sudo apt install gnome-session gdm3
   13  sudo apt --fix-broken install
   14  sudo apt install gnome-session gdm3
   15  gedit &
   16  sudo apt install gedit
   17  gedit &
   18  cd proj
   19  ls
   20  cd hermes
   21  ls
   22  autoconf
   23  sudo apt install autoconf
   24  autoconf
   25  sudo apt install perl
   26  autoconf
   27  ls
   28  cd proj/solar
   29  git pull
   30  cd
   31  sudo apt update
   32  sudo apt install -y ubuntu-desktop xrdp
   33  sudo vi /etc/xrdp/xrdp.ini
   34  sudo gedit /etc/xrdp/xrdp.ini &
   35  sudo service xrdp restart
   36  desktop
   37  rdesktop
   38  history
   39  sudo systemctl enable xrdp









Setting up ui based ubuntu machine on AWS.
In security group open port 5901. Then ssh to the server instance. Run following commands to install ui and vnc server:

sudo apt-get update
sudo apt-get install ubuntu-desktop
sudo apt-get install vnc4server
Then run following commands and enter the login password for vnc connection:

su - awsgui

vncserver

vncserver -kill :1

vim /home/awsgui/.vnc/xstartup
Then hit the Insert key, scroll around the text file with the keyboard arrows, and delete the pound (#) sign from the beginning of the two lines under the line that says "Uncomment the following two lines for normal desktop." And on the second line add "sh" so the line reads

exec sh /etc/X11/xinit/xinitrc. 
When you're done, hit Ctrl + C on the keyboard, type :wq and hit Enter.

Then start vnc server again.

vncserver
You can download xtightvncviewer to view desktop(for Ubutnu) from here https://help.ubuntu.com/community/VNC/Clients

In the vnc client, give public DNS plus ":1" (e.g. www.example.com:1). Enter the vnc login password. Make sure to use a normal connection. Don't use the key files.

Additional guide available here: http://www.serverwatch.com/server-tutorials/setting-up-vnc-on-ubuntu-in-the-amazon-ec2-Page-3.html

Mac VNC client can be downloaded from here: https://www.realvnc.com/en/connect/download/viewer/

Port opening on console

sudo iptables -A INPUT -p tcp --dport 5901 -j ACCEPT

If the grey window issue comes. Mostly because of ".vnc/xstartup" file on different user. So run the vnc server also on same user instead of "awsgui" user.

vncserver

shareimprove this answer
edited Dec 13 '18 at 1:47
answered Sep 5 '14 at 18:32

sugunan
3,54152956
13
I think you're missing su - awsgui after sudo usermod -aG admin awsgui – Konstantin K Dec 11 '14 at 10:38
12
Remember to open port 5901 in your Security Group for this to work. Thanks for the detailed answer! – Daniel Magliola Jan 20 '15 at 11:46
11
Tried a couple of guides, including this one, and I only get a grey background - no ubuntu desktop. – Wrench Feb 11 '15 at 23:10
5
I did this which solved the grey background for me digitalocean.com/community/questions/… – timhc22 Feb 23 '15 at 0:57
7
Try vim .vnc/xstartup if vim awsgui/.vnc/xstartup didn't work – TomasVeras Oct 29 '15 at 0:12
show 13 more comments
 
66

So I follow first answer, but my vnc viewer gives me grey screen when I connect to it. And I found this Ask Ubuntu link to solve that.

The only difference with previous answer is you need to install these extra packages:

apt-get install gnome-panel gnome-settings-daemon metacity nautilus gnome-terminal
And use this ~/.vnc/xstartup file:

#!/bin/sh

export XKL_XMODMAP_DISABLE=1
unset SESSION_MANAGER
unset DBUS_SESSION_BUS_ADDRESS

[ -x /etc/vnc/xstartup ] && exec /etc/vnc/xstartup
[ -r $HOME/.Xresources ] && xrdb $HOME/.Xresources
xsetroot -solid grey
vncconfig -iconic &

gnome-panel &
gnome-settings-daemon &
metacity &
nautilus &
gnome-terminal &
Everything else is the same.

Tested on EC2 Ubuntu 14.04 LTS.

shareimprove this answer
edited Apr 13 '17 at 12:22

Community♦
11
answered Mar 21 '16 at 3:07

yuchien
1,0121913
2
This worked for me on top of the previous answer and with su - awsgui done before running the vnc commands. – Vincenzo Pii Apr 15 '16 at 12:44
1
This step was necessary to get it working after following most tutorials about how to setup ubuntu desktop on aws with tightvncserver. None of the tutorials worked for me without this step. – techdog Jul 15 '16 at 2:37
1
You may need to reboot your OS after following these steps.. I followed this answer and has to reboot first. – tno2007 Nov 24 '16 at 11:26
I tried the above steps and I can see the Ubuntu on RealVNC. But I can see only the terminal and desktop. Somehow other UI parts like Toolbar, Applications etc are missing. Any other steps are there or any fix for this? – Vinayak Jan 22 '17 at 10:47
Run the following in terminal: killall gnome-panel && sudo gnome-panel & – Octocat Apr 29 '17 at 18:12 
show 2 more comments
 
14

For Ubuntu 16.04
1) Install packages

$ sudo apt update;sudo apt install --no-install-recommends ubuntu-desktop
$ sudo apt install gnome-panel gnome-settings-daemon metacity nautilus gnome-terminal vnc4server
2) Edit /usr/bin/vncserver file and modify as below

Find this line

"# exec /etc/X11/xinit/xinitrc\n\n".
And add these lines below.

"gnome-session &\n".
"gnome-panel &\n".
"gnome-settings-daemon &\n".
"metacity &\n".
"nautilus &\n".
"gnome-terminal &\n".
3) Create VNC password and vnc session for the user using "vncserver" command.

lonely@ubuntu:~$ vncserver
You will require a password to access your desktops.
Password:
Verify:
xauth: file /home/lonely/.Xauthority does not exist
New 'ubuntu:1 (lonely)' desktop is ubuntu:1
Creating default startup script /home/lonely/.vnc/xstartup
Starting applications specified in /home/lonely/.vnc/xstartup
Log file is /home/lonely/.vnc/ubuntu:1.log
Now you can access GUI using IP/Domain and port 1

stackoverflow.com:1

Tested on AWS and digital ocean .

For AWS, you have to allow port 5901 on firewall

To kill session

$ vncserver -kill :1
Refer:

https://linode.com/docs/applications/remote-desktop/install-vnc-on-ubuntu-16-04/

Refer this guide to create permanent sessions as service

http://www.krizna.com/ubuntu/enable-remote-desktop-ubuntu-16-04-vnc/

shareimprove this answer







