# @Copyright@
# 
# 				Rocks(r)
# 		         www.rocksclusters.org
# 		            version 5.0 (V)
# 
# Copyright (c) 2000 - 2008 The Regents of the University of California.
# All rights reserved.	
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright
# notice unmodified and in its entirety, this list of conditions and the
# following disclaimer in the documentation and/or other materials provided 
# with the distribution.
# 
# 3. All advertising and press materials, printed or electronic, mentioning
# features or use of this software must display the following acknowledgement: 
# 
# 	"This product includes software developed by the Rocks(r)
# 	Cluster Group at the San Diego Supercomputer Center at the
# 	University of California, San Diego and its contributors."
# 
# 4. Except as permitted for the purposes of acknowledgment in paragraph 3,
# neither the name or logo of this software nor the names of its
# authors may be used to endorse or promote products derived from this
# software without specific prior written permission.  The name of the
# software includes the following terms, and any derivatives thereof:
# "Rocks", "Rocks Clusters", and "Avalanche Installer".  For licensing of 
# the associated name, interested parties should contact Technology 
# Transfer & Intellectual Property Services, University of California, 
# San Diego, 9500 Gilman Drive, Mail Code 0910, La Jolla, CA 92093-0910, 
# Ph: (858) 534-5815, FAX: (858) 534-7345, E-MAIL:invent@ucsd.edu
# 
# THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS''
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS
# BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
# BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
# OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
# IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# @Copyright@

PKGROOT = /opt/mgltools
TARNAME = mgltools_x86_64Linux2_nightly_20110111
DIRNAME = mgltools_x86_64Linux2_latest
MGLBASEPATH = $(PKGROOT)/MGLToolsPckgs

REDHAT.ROOT = $(PWD)/../../
-include $(ROCKSROOT)/etc/Rules.mk
include Rules.mk


install::
	mkdir -p $(ROOT)/$(PKGROOT)
	mkdir -p $(PKGROOT)
	mount --bind $(ROOT)/$(PKGROOT) $(PKGROOT);
	( cd $(DIRNAME); ./install.sh -d $(PKGROOT) -c 1 )
	#remove libstdc from mgltools rpm package if not
	#yum is going crazy and is not installing the system 
	#stdc++
	rm $(PKGROOT)/lib/libstdc++.so.6*
	#apply patch to workflow
	rm -rf $(MGLBASEPATH)/CADD/workflows
	rm -rf $(MGLBASEPATH)/AutoDockTools/VisionInterface/Adt/Macro
	cd patch-files && find . -type f | grep -v CVS | cpio -pduv $(MGLBASEPATH)/CADD
	#script to change server used by opal node in vision
	cp renameOpalService.sh $(PKGROOT)/bin/renameOpalService.sh
	chmod +x $(PKGROOT)/bin/renameOpalService.sh
	umount $(PKGROOT)

build: 
	tar zxvf $(TARNAME).tar.gz

clean::
	rm -rf $(NAMEVERSION)
