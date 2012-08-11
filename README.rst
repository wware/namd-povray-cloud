NAMD and POV-Ray in the Cloud
=============================

This is an attempt to make it easy and painless to do large-scale molecular
modeling and animation using POV-Ray and the NAMD modeling software on VMs in
the cloud, either Google Compute Engine or Amazon Elastic Cloud (EC2). Since
both offer Ubuntu VMs, the goal is to create a debian package and a few shell
scripts which can be pushed to a batch of Ubuntu instances, and then it should
be easy to send jobs to the thing from your desktop and get back results.

Google Compute Engine info
--------------------------

* `Compute Engine product homepage`_
* `Compute Engine developer info`_
* `GCE announcement at IO 2012 keynote day 2`_
* `GCE intro at IO 2012`_
* `GCE details at IO 2012`_
* `Marc Cohen's slides`_
* Comparisons with Amazon offering

  - `Pricing`_
  - `Video transcoding`_

.. _`Compute Engine product homepage`: http://cloud.google.com/products/compute-engine.html
.. _`Compute Engine developer info`: https://developers.google.com/compute/
.. _`GCE announcement at IO 2012 keynote day 2`: http://www.youtube.com/watch?v=tPtJd6AzU8c#t=39m
.. _`GCE intro at IO 2012`: http://www.youtube.com/watch?v=0-sF5ZWB_FY
.. _`GCE details at IO 2012`: http://www.youtube.com/watch?v=ws2VRHq5ars
.. _`Marc Cohen's slides`: http://commondatastorage.googleapis.com/marc-pres/gce-0812/index.html
.. _`Pricing`: http://blog.abourget.net/2012/6/28/amazon-elastic-cloud-computing-ec2-google-cloud-compute-engine/
.. _`Video transcoding`: http://video.heidisoft.com/blog/first-look-google-compute-engine-video-transcoding

Amazon Elastic Cloud info
-------------------------

* `Elastic Cloud product homepage`_
* `EC2 Getting Started Guide`_
* `EC2 User Guide`_
* `EC2 API Reference`_
* `EC2 Technical FAQ`_
* `AWS Marketplace`_ (deploy-ready EC2 apps for sale)

.. _`Elastic Cloud product homepage`: http://aws.amazon.com/ec2/
.. _`EC2 Getting Started Guide`: http://docs.amazonwebservices.com/AWSEC2/latest/GettingStartedGuide/Welcome.html
.. _`EC2 User Guide`: http://docs.amazonwebservices.com/AWSEC2/latest/UserGuide/Welcome.html
.. _`EC2 API Reference`: http://docs.amazonwebservices.com/AWSEC2/latest/APIReference/Welcome.html
.. _`EC2 Technical FAQ`: http://docs.amazonwebservices.com/AWSEC2/latest/UserGuide/TechnicalFAQ.html
.. _`AWS Marketplace`: https://aws.amazon.com/marketplace/ref=mkt_ste_ec2
