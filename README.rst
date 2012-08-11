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
