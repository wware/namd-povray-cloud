cp d45.tgz /usr/local/d45.tgz
cd /usr/local
tar xzvf d45.tgz
rm -f /usr/bin/pyvon
ln -s /usr/local/pyvon/pyvon /usr/bin/pyvon
rm -f /usr/local/d45.tgz
