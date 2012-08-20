Prototyping this setup on VBox
==============================

I want to set up the thing with two disk partitions on VBox. To do that, you
need to use a VMDK disk for the second disk, which will then appear as
/dev/sdb. Then::

 fdisk /dev/sdb
 # Set up one partition for the whole VMDK disk, of type 83 (Linux)
 sudo mkdir /db
 sudo chmod 777 /db
 sudo mkfs.ext4 /dev/sdb1
 sudo mount -t ext4 /dev/sdb1 /db

To change the data directory for MySQL, edit /etc/mysql/my.cnf [1_] and also
tweak /etc/init.d/apparmor and some directory permissions and ownerships [2_].
A similar hack [3_, 4_] for MongoDB.

.. _1: http://niblets.wordpress.com/2007/03/23/changing-your-mysql-data-directory/
.. _2: http://ubuntuforums.org/showpost.php?p=5220182&postcount=2
.. _3: http://www.mongodb.org/display/DOCS/A+Sample+Configuration+Session
.. _4: http://stackoverflow.com/questions/5961145/changing-mongodb-data-storing-directory

It would be good to have any of these post-install mods done in a script, and
eventually put that script into a debian package.
