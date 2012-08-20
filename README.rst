NAMD and POV-Ray in the Cloud
=============================

This is an attempt to make it easy and painless to do large-scale molecular
modeling and animation using POV-Ray and the NAMD modeling software on VMs in
the cloud, either Google Compute Engine or Amazon Elastic Cloud (EC2). Since
both offer Ubuntu VMs, the goal is to create a debian package and a few shell
scripts which can be pushed to a batch of Ubuntu instances, and then it should
be easy to send jobs to the thing from your desktop and get back results.

I'm starting with EC2 because it's available immediately. I'll figure out GCE
later.

Both use Ubuntu 12.04, so I'll need to upgrade my own efforts from Lucid which
I've used in the past.

Hand-wringing
-------------

From [1_]:

 Virtual cores are useless for GROMACS because it saturates the processor
 pipelines on its own... You need real cores, and to configure MPI to assign
 processes 1-to-1 to real cores. The details will be up to you and your system
 admins. You don't need to do anything with GROMACS other than configure with
 MPI.

Gigabit ethernet routers are cheap. Think about the earlier plan of picking up
some quad-core machines, like this [2_]. A couple of them would cost $700 and
would give me eight cores connected by gigabit ethernet. Each has a 500 GB
hard drive and 4 GB of RAM. That should be OK if I don't run X Windows on
them. But the big HD is suggestive of a database cluster, so this might be the
setup for the autosci stuff.

.. _1: http://comments.gmane.org/gmane.science.biology.gromacs.user/44432
.. _2: http://www.provantage.com/acer-pt-sg9p2-003~7ACED0J9.htm

How to make this thing accessible remotely? I won't be running a
server out of my house, but I can use the eApps machine as a
reflector. I guess a web API is the way to go.

Envision four of these machines as the Clojure cluster I was thinking
about earlier, but more thoughtfully optimized.

Let's partition each hard drive into two pieces. One is for booting,
and can easily be overwritten with a new Ubuntu ISO.

The four non-booting partitions could be used as a RAID5 array [3_, 4_]. It
looks like the latency and general performance of RAID5 is pretty bad, but
maybe it's acceptable if you're doing large-ish block writes and reads, and I
bet that's why Google does things 64 MBs at a time. It looks like Ubuntu's
software RAID5 only works if all the drives are run by one motherboard, so
it's inapplicable.

.. _3: http://en.wikipedia.org/wiki/Standard_RAID_levels#RAID_5
.. _4: https://help.ubuntu.com/community/Installation/SoftwareRAID

Maybe HDFS [5_] would be a better choice. Since Hadoop is written in Java
anyway, Hadoop/HDFS is a good fit with the Clojure idea.

.. _5: http://en.wikipedia.org/wiki/Apache_Hadoop#Hadoop_Distributed_File_System

Maybe Gromacs on EC2 is not necessarily so bad. A four-core machine can be
thought of as four CPUs sharing a rack. Two or three of those, and you have
most of the internal communications on the same die, which is better than
gigabit ethernet.

It would cost hundreds of dollars to build this cluster. If it costs only tens
of dollars to make a simulation/animation on EC2, and I only do a few, that's
a win. It would be great to figure out how to make some money off these
animations.

There are plenty of working distributed databases and any of them could be
used in essentially the same way for parallel computing. Any of these would
work: MySQL, PostgreSQL, Hadoop/HBase, MongoDB, most NoSQL [6_] databases.

.. _6: http://en.wikipedia.org/wiki/NoSQL.

It would be very easy to write a tuple-space-like distributed computing system
on MySQL Cluster [7_]. But don't store large blobs in MySQL, store them in
something else. I think it makes sense to have a combination of different
databases: MySQL Cluster, HBase, and MongoDB.

.. _7: http://en.wikipedia.org/wiki/MySQL_Cluster

Assuming I go with my own hardware, the two-partition scheme allows me to
change the distro running on the machines while retaining all the DB content.
Can I do the same with VMs running on a single machine, or with instances in
a cloud? I think both are doable. That would be a great way to set the thing
up cheaply and get it started, and then migrate to EC2 when it becomes large
and unwieldy. I'll experiment with Ubuntu 12.04 server on VBox.
