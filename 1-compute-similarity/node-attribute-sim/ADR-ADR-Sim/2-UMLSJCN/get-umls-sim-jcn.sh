# Bash Shell Script for Computing UMLS-based ADR-ADR Similarity (JCN)
#
# The following shell script requires the UMLS database (MySQL)
# and the Perl module UMLS::Similarity installed first.
#
# This task was completed using the high-performance computing server of
# CBDD Group, Central South University. It took about 48 hours to complete.
#
# Author: Nan Xiao <me@nanx.me>
#
# Date: Aug 19, 2013

umls-similarity.pl --infile ADR.txt --matrix --config queryconfig --measure jcn --username root --password your_password --socket /var/lib/mysql/mysql.sock --database yourumlsdatabase >> umls-sim-jcn.txt
