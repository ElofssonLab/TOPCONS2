<?xml version="1.0"?>
<!DOCTYPE book PUBLIC "-//OASIS//DTD DocBook XML V4.5//EN"
                    "http://www.oasis-open.org/docbook/xml/4.5/docbookx.dtd"
[
<!ENTITY % MATHML.prefixed "INCLUDE">
<!ENTITY % MATHML.prefix "mml">
<!ENTITY % equation.content "(alt?, (graphic+|mediaobject+|mml:math))">
<!ENTITY % inlineequation.content 
                "(alt?, (inlinegraphic+|inlinemediaobject+|mml:math))">
<!ENTITY % mathml PUBLIC "-//W3C//DTD MathML 2.0//EN"
        "http://www.w3.org/Math/DTD/mathml2/mathml2.dtd">
]>


<article xmlns='http://docbook.org/ns/docbook'
         xmlns:mml="http://www.w3.org/1998/Math/MathML">
  <articleinfo>
    <title>modhmm</title>
  </articleinfo>
  <sect1 id="introduction">
    <title>Introduction</title>
    <para>
<mediaobject>
  <imageobject>
    <imagedata fileref="images/modhmm_web_pic.jpg" format="JPEG"/>
  </imageobject>
</mediaobject>



modhmm is software tool for building, training and scoring hidden Markov models. The software is open source and licensed under the GPL license.
Building a hidden Markov model is done using a modular approach which means that there is a set of
predefined model subparts (modules) to choose from when putting together a complete hmm. The idea is to
simplify the building of large hmms while still allowing for hmms of arbitrary architecture to be built.
The model building tool is called modhmmc.
</para><para>
modhmm consists of three main subparts, modhmmc, modhmmt and modhmms for creating training and scoring hmms respectively.
Training and scoring can be done using either single sequences, multiple sequence alignments or sequence profiles as training/scoring input.
</para><para>
The current state of modhmm is still somewhat preliminary.

  </para>
    <para>
        The primary URL for this document is
         <ulink url="http://modhmm.sourceforge.net">http://modhmm.sourceforge.net</ulink>.
    </para>
  </sect1>


  <sect1 id="hmg_file_format">
    <title>The .hmg file format</title>
    <para>
There are two different formats for storing modhmm-models. These formats are in most
aspects identical. Their differences are associated to the use of multiple alphabets or not.
The multiple alphabet format (up to four different parallel alphabets are possible at the time) includes
some additional entries for this information which are not present in the single alphabet format.
Every model created by modhmmc or modhmmt is saved in the <emphasis>.hmg</emphasis> text format.
It is fairly sensitive to minor changes as adding extra blank lines, extra blanks within a line, etc.
This type of changes may work, but nothing is guaranteed. Therefore caution is decreed when manually editing a .hmg file
    </para>

  <sect2 id="the_single_alphabet_format">
    <title>The single alphabet format</title>
  <sect3 id="the_single_alphabet_format_header">
    <title>Header</title>
    <para>
The header of a .hmg file contains 14 lines (+ 2 compulsory blank lines).
<programlisting><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/example_files/hmg-header">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include></programlisting>

<variablelist><title>Description of headers</title>
<varlistentry><term>NAME</term><listitem><para>the name of the model
</para></listitem></varlistentry><varlistentry><term>TIME OF CREATION</term><listitem><para>time for last modification of the model
</para></listitem></varlistentry><varlistentry><term>ALPHABET</term><listitem><para>the alphabet of the emissions, each letter separated by ';'
</para></listitem></varlistentry><varlistentry><term>ALPHABET LENGTH</term><listitem><para>the number of letters in the alphabet
</para></listitem></varlistentry><varlistentry><term>NR OF MODULES</term><listitem><para>the number of hmm modules
</para></listitem></varlistentry><varlistentry><term>NR OF VERTICES</term><listitem><para>the number of states
</para></listitem></varlistentry><varlistentry><term>NR OF TRANSITIONS</term><listitem><para>the total number of transitions
</para></listitem></varlistentry><varlistentry><term>NR OF DISTRIBUTION GROUPS</term><listitem><para>the number of distribution groups. A distribution group is a set of states whose emission
  probabilities have been tied together so that during training updating, they are regarded as the same state.
</para></listitem></varlistentry><varlistentry><term>NR OF TRANSITION TIE GROUPS</term><listitem><para>the number of tied transitions. A transition tie is the same thing as a distribution group, but for
  transitions.
</para></listitem></varlistentry><varlistentry><term>NR OF EMISSION PRIORFILES</term><listitem><para>number of emission priorfiles
</para></listitem></varlistentry><varlistentry><term>EMISSION PRIORFILES</term><listitem><para>names (and paths) of emission priorfiles. An emission priorfile is a file with prior information over
  the emissions used to weight the observed emission frequences against a belief prior to the observation.
</para></listitem></varlistentry><varlistentry><term>NR OF TRANSITION PRIORFILES</term><listitem><para>number of transition priorfiles
</para></listitem></varlistentry><varlistentry><term>TRANSITION PRIORFILES</term><listitem><para>names (and paths) of transition priorfiles. Same thing as an emission priorfile, but for transitions. 
</para></listitem></varlistentry>
</variablelist>


    </para>

  </sect3>
  <sect3 id="the_single_alphabet_format_modules">
    <title>Modules</title>
    <para>
Each module section contains 4 rows in the beginning (+ a compulsory blank row) and a set of vertex sections, each separated
by a blank row. Each module section is ended by a blank row and a row of '-'.

<programlisting><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/example_files/hmg-modules">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include></programlisting>


<variablelist><title>Description of modules variables</title>
<varlistentry><term>Module</term><listitem><para>the name of the module
</para></listitem></varlistentry><varlistentry><term>Type</term><listitem><para>the type of the module ( see <xref linkend="modules"/> )
</para></listitem></varlistentry><varlistentry><term>NrVertices</term><listitem><para>number of states in this module
</para></listitem></varlistentry><varlistentry><term>Emission prior file</term><listitem><para>possible emission prior file associated with this module, 'null' means that no file is associated 
</para></listitem></varlistentry><varlistentry><term>Transition prior file</term><listitem><para>possible transition prior file associated with this module, 'null' means that no file is associated 
</para></listitem></varlistentry></variablelist>
    </para>
  </sect3>

  <sect3 id="the_single_alphabet_format_vertices">
    <title>Vertices</title>
    <para>

Each vertex section consists of 9 initial rows, followed by a sections for transition probabilities, end transition probabilities
and emission probabilities respectively.
    </para>


<programlisting>Vertex 1:
Vertex type: standard
Vertex label: d
Transition prior scaler: 1.0
Emission prior scaler: 1.0
Nr transitions = 1
Nr end transitions = 0
Nr emissions = 20
Transition probabilities
        Vertex 2: 1.0
End transition probabilities
Emission probabilities
        A: 0.05
        C: 0.05
        .
        .
        .
        T: 0.05
        V: 0.05
        W: 0.05
        Y: 0.05
</programlisting>



<variablelist 
><varlistentry><term 
>
Vertex </term><listitem>
     <!--l. 115--><para>the number of the state
     </para></listitem></varlistentry><varlistentry><term 
>
Vertex type </term><listitem>
     <!--l. 116--><para>the type of the state. 3 types exist: standard, silent and locked. A
     standard state is an regular emitting state. A silent state is a state that
     does not emit any symbols. Finally a locked state is an emitting state
     for which the emission probabilities are fixed, i.e. will not be updated
     during training.
     </para></listitem></varlistentry><varlistentry><term 
>
Transition prior scaler </term><listitem>
     <!--l. 119--><para>the transition prior scaler is a factor which describes how much weight
     to put on the possible prior distribution associated with this vertex

     </para></listitem></varlistentry><varlistentry><term 
>
Emission prior scaler </term><listitem>
     <!--l. 121--><para>the emission prior scaler is a factor which describes how much weight to
     put on the possible Dirichlet prior mixture associated with this vertex.
     </para></listitem></varlistentry><varlistentry><term 
>
Nr transitions </term><listitem>
     <!--l. 123--><para>is the number of transitions from this state (except to end states)
     </para></listitem></varlistentry><varlistentry><term 
>
Nr end transitions </term><listitem>
     <!--l. 124--><para>is the number of transitions from this state to end states
     </para></listitem></varlistentry><varlistentry><term 
>
Nr emissions </term><listitem>
     <!--l. 125--><para>is the number of different possible emissions from this state, usually
     equal to the alphabet size.</para></listitem></varlistentry></variablelist>


The next three are all listings of transition, end transition and emission
probabilities respectively.
     <variablelist 
><varlistentry><term 
>
Transition probabilities </term><listitem>
     <!--l. 129--><para>has  a  row  for  each  of  the  transitions  which  states  which  state  the
     transition is to and the probability of the transition</para></listitem></varlistentry></variablelist>

<programlisting>Transition probabilities
        Vertex 2: 0.4
        Vertex 3: 0.6
</programlisting>


     <variablelist 
><varlistentry><term>End transition probabilities</term><listitem><para>
 same as above, but for end transitions, usually there are not more than one of this type of
  transition, but nothing in the program prohibits this   </para></listitem></varlistentry></variablelist>


<programlisting>End transition probabilities
        Vertex 4: 1.0
</programlisting>


     <variablelist 
><varlistentry><term>Emission probabilities</term><listitem><para>  the probabilities for emitting the different letters in the alphabet,
  one row for each letter in the alphabet. In the case of continuous emissions, the alphabet is interpreted
  by the program as follows. The letters are divided into groups of 3. Each group describe a mixture component
  of a mixture of one-dimensional normal distributions.
  For each group, the first letter represents the mean value, the second the variance and the third the coefficient
  for the particular mixture component.
     </para></listitem></varlistentry></variablelist>


Discrete:
<programlisting>Emission probabilities
        A: 0.05
        C: 0.05
        D: 0.05
        E: 0.05
        .
        .
        .
        Y: 0.05
</programlisting>

Continuous (2 mixture components):
<programlisting>Emission probabilities
        m1: 22.5
        var1: 0.34
        co1: 0.45
        m2: -22.3
        var2: 3.11
        co2: 0.55
</programlisting>






  </sect3>
  <sect3 id="the_single_alphabet_format_emission_distribution_groups">
    <title>Emission distribution groups</title>
    <para>
The emission distribution group  has a line for each distribution group which simply states the numbers of the vertices
that belong to a particular group

<programlisting>Group 1: 1 2
Group 2: 3 4 
</programlisting>

    </para>
  </sect3>
  <sect3 id="the_single_alphabet_format_transition_tie_groups">
    <title>Transition tie groups</title>
    <para>
The transition tie groups entries (one line) describe two or more transitions which are tied together. A transition is specified
as <emphasis>from-state arrow to-state</emphasis>.

<programlisting>Tie 1: 1->2 3->4
Tie 2: 5->6 7->10000 8->10001
</programlisting>



    </para>
  </sect3>




  </sect2>
  </sect1>












  <sect1 id="software">
    <title>Software</title>


  <sect2 id="download">
    <title>Download</title>
    <para>      
       Download the software from the <ulink url="http://sourceforge.net/projects/modhmm">sourceforge</ulink> project page. 

 The latest version of modhmm is @PACKAGE_VERSION@. 
    </para>
  </sect2>
  <sect2 id="installation">
    <title>Installation</title>

    <sect3 id="installation_with_prebuilt_package">
      <title>Installation with prebuilt package</title>

    <sect4 id="installation_on_windows">
      <title>Installation on Windows with .exe file</title>
<para>
To install modhmm on Windows, first download the modhmm-@PACKAGE_VERSION@-win32.exe and then execute the file ( click on it ). 
    </para>
  </sect4>


    <sect4 id="installation_on_ubuntu_and_debian">
      <title>Installation on Ubuntu and Debian</title>

<para>
To install modhmm on Ubuntu or Debian, first download the modhmm-@PACKAGE_VERSION@.deb  and then log in as root and  
<programlisting><![CDATA[
# dpkg -i modhmm-@PACKAGE_VERSION@.deb 
]]></programlisting>

    </para>
  </sect4>


    <sect4 id="installation_on_centos_and_fedora_linux">
      <title>Installation on CentOS and Fedora</title>

<para>
To install modhmm on Centos or Debian, first download the modhmm-@PACKAGE_VERSION@.Linux.rpm   and then log in as root and  
<programlisting><![CDATA[
# yum localinstall modhmm-@PACKAGE_VERSION@.Linux.rpm 
]]></programlisting>

    </para>
  </sect4> 





    <sect4 id="installation_MacOS_X">
      <title>Installation on Mac OS X</title>

<para>
To install modhmm on a Mac OS X v10.5 ( Leopard ) on a Mac computer with Intel cpu, first download the modhmm-@PACKAGE_VERSION@-MacOSX10.5.tar.gz   and then 
<programlisting><![CDATA[
$ tar xfz modhmm-@PACKAGE_VERSION@-MacOSX10.5.tar.gz  
]]></programlisting>

    </para>
<para>
To install modhmm on a Mac OS X v10.4 ( Tiger ) on a Mac computer with Intel cpu, first download the modhmm-@PACKAGE_VERSION@-MacOSX10.4.tar.gz   and then 
<programlisting><![CDATA[
$ tar xfz modhmm-@PACKAGE_VERSION@-MacOSX10.4.tar.gz
]]></programlisting>

    </para>
  </sect4>




 </sect3>




    <sect3 id="building_from_source">
      <title>Building from source</title>

    <sect4 id="building_from_source_on_unix">
      <title>Building from source on Unix</title>

      <para>To build modhmm on Unix ( e.g. Linux, MacOSX, CygWin ) you need to have this installed
        <itemizedlist mark="bullet">
          <listitem>
            <para>
              <ulink url="http://www.cmake.org">cmake</ulink>
            </para>
          </listitem>
          <listitem>
            <para>
              <ulink url="http://xmlsoft.org/">libxml2</ulink>
            </para>
          </listitem>
        </itemizedlist>
      </para>
      <para>

If you have the modhmm source code in the directory <filename>/tmp/modhmm</filename> and you want to install modhmm into the directory <filename>/tmp/install</filename>, you 

First run <command>cmake</command> then <command>make</command> and then <command>make install</command>
<programlisting><![CDATA[
$ mkdir /tmp/build
$ cd /tmp/build
$ cmake -DCMAKE_INSTALL_PREFIX=/tmp/install /tmp/source && make && make install
-- Configuring done
-- Generating done
-- Build files have been written to: /tmp/build
Scanning dependencies of target fastdist
[  3%] Building CXX object src/c++/CMakeFiles/fastdist.dir/programs/fastdist.o
[  6%] Building CXX object src/c++/CMakeFiles/fastdist.dir/BitVector.o
[  9%] Building CXX object src/c++/CMakeFiles/fastdist.dir/Exception.o
[ 12%] Building CXX object src/c++/CMakeFiles/fastdist.dir/InitAndPrintOn_utils.o
[ 15%] Building CXX object src/c++/CMakeFiles/fastdist.dir/Object.o
[ 18%] Building CXX object src/c++/CMakeFiles/fastdist.dir/Sequence.o
[ 21%] Building CXX object src/c++/CMakeFiles/fastdist.dir/SequenceTree.o
[ 25%] Building CXX object src/c++/CMakeFiles/fastdist.dir/SequenceTree_MostParsimonious.o
[ 28%] Building CXX object src/c++/CMakeFiles/fastdist.dir/Simulator.o
[ 31%] Building CXX object src/c++/CMakeFiles/fastdist.dir/arg_utils_ext.o
[ 34%] Building CXX object src/c++/CMakeFiles/fastdist.dir/file_utils.o
[ 37%] Building CXX object src/c++/CMakeFiles/fastdist.dir/stl_utils.o
[ 40%] Building CXX object src/c++/CMakeFiles/fastdist.dir/DNA_b128/DNA_b128_String.o
[ 43%] Building CXX object src/c++/CMakeFiles/fastdist.dir/DNA_b128/Sequences2DistanceMatrix.o
[ 46%] Building CXX object src/c++/CMakeFiles/fastdist.dir/aml/AML_LeafLifting.o
[ 50%] Building CXX object src/c++/CMakeFiles/fastdist.dir/aml/AML_given_edge_probabilities.o
[ 53%] Building CXX object src/c++/CMakeFiles/fastdist.dir/aml/AML_local_improve.o
[ 56%] Building CXX object src/c++/CMakeFiles/fastdist.dir/aml/AML_star.o
[ 59%] Building CXX object src/c++/CMakeFiles/fastdist.dir/aml/Big_AML.o
[ 62%] Building CXX object src/c++/CMakeFiles/fastdist.dir/distance_methods/LeastSquaresFit.o
[ 65%] Building CXX object src/c++/CMakeFiles/fastdist.dir/distance_methods/NeighborJoining.o
[ 68%] Building CXX object src/c++/CMakeFiles/fastdist.dir/sequence_likelihood/Kimura2parameter.o
[ 71%] Building CXX object src/c++/CMakeFiles/fastdist.dir/sequence_likelihood/TamuraNei.o
[ 75%] Building CXX object src/c++/CMakeFiles/fastdist.dir/sequence_likelihood/ambiguity_nucleotide.o
[ 78%] Building CXX object src/c++/CMakeFiles/fastdist.dir/sequence_likelihood/dna_pairwise_sequence_likelihood.o
[ 81%] Building CXX object src/c++/CMakeFiles/fastdist.dir/sequence_likelihood/string_compare.o
[ 84%] Building CXX object src/c++/CMakeFiles/fastdist.dir/DistanceMatrix.o
[ 87%] Building C object src/c++/CMakeFiles/fastdist.dir/arg_utils.o
cc1: warning: command line option "-fno-default-inline" is valid for C++/ObjC++ but not for C
[ 90%] Building C object src/c++/CMakeFiles/fastdist.dir/std_c_utils.o
cc1: warning: command line option "-fno-default-inline" is valid for C++/ObjC++ but not for C
[ 93%] Building C object src/c++/CMakeFiles/fastdist.dir/DNA_b128/sse2_wrapper.o
[ 96%] Building CXX object src/c++/CMakeFiles/fastdist.dir/DNA_b128/computeTAMURANEIDistance_DNA_b128_String.o
[100%] Building CXX object src/c++/CMakeFiles/fastdist.dir/DNA_b128/computeDistance_DNA_b128_String.o
Linking CXX executable fastdist
[100%] Built target fastdist
[100%] Built target fastdist
Linking CXX executable CMakeFiles/CMakeRelink.dir/fastdist
Install the project...
-- Install configuration: ""
-- Install configuration: ""
-- Installing /tmp/install/bin/fastdist
-- Install configuration: ""
]]></programlisting>

If you want to build the html documentation ( i.e. this page ) you need to pass the -DBUILD_DOCBOOK=ON option to <application>cmake</application>.
      </para>
    </sect4>


    <sect4 id="building_install_packages">
      <title>Building install packages</title>
      <para>This is section is mainly intended for package maintainers</para>

    <sect5 id="building_an_exe_file_for_windows">
      <title>Building an .exe file for Windows</title>
<para>
To build the modhmm nullsoft installer package ( modhmm-@PACKAGE_VERSION@-win32.exe ) 
you need to have this installed
        <itemizedlist mark="bullet">
          <listitem>
            <para>
              <ulink url="http://www.cmake.org">cmake</ulink>
            </para>
          </listitem>
          <listitem>
            <para>
              <ulink url="http://www.mingw.org">mingw</ulink>
            </para>
          </listitem>
          <listitem>
            <para>
              <ulink url="http://www.mingw.org/wiki/msys">msys</ulink>
            </para>
          </listitem>

          <listitem>
            <para>
              <ulink url="http://gnuwin32.sourceforge.net/packages/wget.htm">wget</ulink>
            </para>
          </listitem>

          <listitem>
            <para>
              <ulink url="http://nsis.sourceforge.net">Nullsoft Scriptable Install System</ulink>
            </para>
          </listitem>
        </itemizedlist>
on your Windows machine.
      </para>
      <para>

Just open up a msys bash shell 
<programlisting><![CDATA[
$ mkdir tmpbuild
$ cd tmpbuild
$ cmake path/to/the/modhmm/source/code  -DSTATIC=ON -G "MSYS Makefiles" && make win32installer
]]></programlisting>

The source code for gengetopt, libz and libxml will be automatically downloaded and built statically.
      </para>
    </sect5>


    <sect5 id="building_install_package_rpm">
      <title>Building an rpm</title>
<para>
On a CentOS or Fedora machine, first log in as root and install the dependencies
<programlisting><![CDATA[
# yum install xmlto libxml2-devel cmake gcc-c++ binutils gengetopt
]]></programlisting>

Check that cmake is version 2.6 or later
<programlisting><![CDATA[
$ cmake --version
cmake version 2.6-patch 0
]]></programlisting>
If it is older you could download a cmake binary directly from <ulink url="http://www.cmake.org">www.cmake.org</ulink>

<programlisting><![CDATA[
$ mkdir /tmp/build
$ cd /tmp/build
$ cmake -DCMAKE_INSTALL_PREFIX=/ -DBUILD_DOCBOOK=ON /tmp/source && make package
]]></programlisting>

      </para>
    </sect5>

    <sect5 id="building_install_package_deb">
      <title>Building a deb package</title>
<para>
On a Debian or Ubuntu machine, first log in as root and install the dependencies
<programlisting><![CDATA[
# apt-get install libxml2-dev cmake g++ binutils gengetopt
]]></programlisting>

Check that cmake is version 2.6 or later
<programlisting><![CDATA[
$ cmake --version
cmake version 2.6-patch 0
]]></programlisting>
If it is older you could download a cmake binary directly from <ulink url="http://www.cmake.org">www.cmake.org</ulink>. Now build the deb package.


<programlisting><![CDATA[
$ mkdir /tmp/build
$ cd /tmp/build
$ cmake -DCMAKE_INSTALL_PREFIX=/ -DBUILD_DOCBOOK=ON /tmp/source && make package
]]></programlisting>

      </para>
    </sect5>


    <sect5 id="building_install_package_for_macosx">
      <title>Building install package for MacOS X</title>
<para>
To build the modhmm install package for MacOS X
you need to have this installed
        <itemizedlist mark="bullet">
          <listitem>
            <para>
              <ulink url="http://www.cmake.org">cmake</ulink>
            </para>
          </listitem>
          <listitem>
            <para>
              <ulink url="http://www.gnu.org/software/wget/">wget</ulink>
            </para>
          </listitem>
        </itemizedlist>
on your MacOS X computer.
      </para>


      <para>

Check that cmake is version 2.6 or later
<programlisting><![CDATA[
$ cmake --version
cmake version 2.6-patch 0
]]></programlisting>


<programlisting><![CDATA[
$ mkdir /tmp/build
$ cd /tmp/build
$ cmake -DSTATIC=ON -DCPACK_GENERATOR="TGZ" /tmp/source && make package
]]></programlisting>

      </para>
    </sect5>




    </sect4>



    </sect3>





  </sect2>



  <sect2 id="usage">
    <title>Usage</title>


  <sect3 id="modhmmc">
    <title>modhmmc</title>

<para>
<application>modhmmc</application> is the program for designing an HMM. The main <application>modhmmc</application> program has a command line based interactive user interface which
lets the user specify alphabet, states, transition and emission probabilities, etc. If an hmm with one alphabet is designed,
it will automatically be saved in the one-alphabet format. If an hmm with multiple alphabets is designed, it will automatically
be saved in the multiple alphabet format
</para>




  <sect4 id="modhmmc_command_line_options">
    <title>Command line options</title>


<para>

Type <userinput>modhmmc --help</userinput> to see the command line options

<programlisting><![CDATA[
[user@saturn ~]$ ]]><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/commands/modhmmc_help">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
<xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_BINARY_DIR}/xincluded_files/modhmmc_help">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
</programlisting>
</para></sect4>





















  <sect4 id="modules">
    <title>Modules</title>
    <para>      
The states of an HMM are specified as a collection of <emphasis>state modules</emphasis> with transitions between them.
A module is a set of states which are interconnected in a predecided fashion. The idea behind modules
is to make the creation of large HMMs easier. HMMs with several hundred states are not uncommon in
sequence analysis and specifying each and every state and transition in such a case is impractical.
modhmmc currently has 7 module types to choose from: singlenode, singleloop, forward std, forward alt, cluster, profile7 and
profile9 ).
modhmmc allows for the creation of both regular and silent states.

</para>

<para>

Transition probabilities are set by default inside the modules to
correspond to the intrinsic properties of that module (see descriptions of the modules). Transition
probabilities between modules are by default set uniformly, so that the probabilities of going from a module to another are
equally distributed among all its neighbours. Emission probabilities for each state are either set manually by the user
or according to a chosen distribution. The choices are to set the values uniformly, randomly or to zero for all
letters, which creates a silent state. There are also some special distributions that are particular to certain states.


</para>


  <sect5 id="singlenode_module">
    <title>Singlenode module</title>


<para>


<mediaobject>
  <imageobject>
    <imagedata fileref="images/singlenodemod.png" format="PNG"/>
  </imageobject>
  <imageobject>
    <imagedata fileref="images/singlenodemod.eps" format="EPS"/>
  </imageobject>
</mediaobject>


The singlenode module is the most basic module of modhmmc. It consists of only one state, emitting by default, but the user
may specify it as silent. All other modules may be built using collections of
singlenodes. The singlenode has no input parameters.
</para><para>
The singlenode module is the most basic module of modhmmc. It consists of only one state, emitting by default, but the user
may specify it as silent. All other modules may be built using collections of
singlenodes. The singlenode has no input parameters.
</para>
</sect5>


  <sect5 id="singleloop_module">
    <title>Singleloop module</title>
<para>



<mediaobject>
  <imageobject>
    <imagedata fileref="images/singleloopmod.png" format="PNG"/>
  </imageobject>
  <imageobject>
    <imagedata fileref="images/singleloopmod.eps" format="EPS"/>
  </imageobject>
</mediaobject>


The singleloop module consists of one state with a transition to itself. It is necessarily emitting, since no loops of
silent states are allowed. The singleloop has one input parameter. The user specifies the expected length of the loop, which
sets the loop transition probability to
<inlineequation>

  <mml:math>
<mml:mrow>
 <mml:msub>
   <mml:mi>p</mml:mi>
   <mml:mn>trans</mml:mn>
 </mml:msub>
 <mml:mo>=</mml:mo>
 <mml:msup>
   <mml:mi>e</mml:mi>
   <mml:mfrac><mml:mrow><mml:mi>ln</mml:mi><mml:mn>0.5</mml:mn></mml:mrow><mml:mi>L</mml:mi></mml:mfrac>
 </mml:msup>
</mml:mrow>
  </mml:math>
  

</inlineequation>


 where 
<inlineequation>
  <mml:math>
<mml:mi>L</mml:mi>
  </mml:math>
</inlineequation>

is the specified loop length. This
results in the probability of a loop length longer than or equal to <inlineequation>
  <mml:math>
<mml:mi>L</mml:mi>
  </mml:math>
</inlineequation>
 being equal to the probability of a loop
length shorter than <inlineequation>
  <mml:math>
<mml:mi>L</mml:mi>
  </mml:math>
</inlineequation>
, that is, <inlineequation>
  <mml:math>
<mml:mi>L</mml:mi>
  </mml:math>
</inlineequation>
 is the expected median loop length. 



</para>
</sect5>


  <sect5 id="forward_std_module">
    <title>Forward std module</title>
<para>



<mediaobject>
  <imageobject>
    <imagedata fileref="images/forwardmod.png" format="PNG"/>
  </imageobject>
  <imageobject>
    <imagedata fileref="images/forwardmod.eps" format="EPS"/>
  </imageobject>
</mediaobject>


The forward module is a set of states connected in a straight line, indexed 1 to <varname>n</varname>. All states are emitting by default.
The two input parameters <varname>(m,n)</varname>
specify the shortest and longest possible routes through the module. The number of states in the module is equal to the
length of the longest possible route <varname>n</varname>. For the state with index <varname>m</varname>-1 there is a transition to all states in the module with
a higher index number up to and including the state with index <varname>n</varname>.
The transition probabilities are by default set so that the total probability is equal for all possible paths through the module,
but it is also possible to set the length probilities according to a binomial distribution, with the base either in the shortest 
or the longest path.
Any state that connects to this module will connect to the state with index 1, and all outgoing connections from this module
go from the state with index <varname>n</varname>.


</para>
</sect5>


  <sect5 id="forward_alt_module">
    <title>Forward alt module</title>
<para>

<mediaobject>
  <imageobject>
    <imagedata fileref="images/forwardmod_alt.png" format="PNG"/>
  </imageobject>
  <imageobject>
    <imagedata fileref="images/forwardmod_alt.eps" format="EPS"/>
  </imageobject>
</mediaobject>


The forward module is a set of states connected in a straight line, indexed 1 to <varname>n</varname>. All states are emitting by default.
The two input parameters <varname>(m,n)</varname>
specify the shortest and longest possible routes through the module. The number of states in the module is equal to the
length of the longest possible route <varname>n</varname>. For all states with index <varname>m</varname>-1  to <varname>n</varname>-2 there is a transition to the last state.
The transition probabilities are set so that the total probability is equal for all possible paths through the module,
but it is also possible to set the length probilities according to a binomial distribution, with the base either in the shortest 
or the longest path.
Any state that connects to this module will connect to the state with index 1, and all outgoing connections from this module
go from the state with index <varname>n</varname>.


</para>
</sect5>


  <sect5 id="cluster_module">
    <title>Cluster module</title>
<para>


<mediaobject>
  <imageobject>
    <imagedata fileref="images/clustermod.png" format="PNG"/>
  </imageobject>
  <imageobject>
    <imagedata fileref="images/clustermod.eps" format="EPS"/>
  </imageobject>
</mediaobject>


The cluster module is a fully interconnected set of states. Every state has a transition to every other state.
The transition probabilities are evenly distributed by default. All states are emitting. The input parameter is the number of states.
Incoming transitions connect to all states, and outgoing connections go from all states.


</para>
</sect5>


  <sect5 id="profile9_module">
    <title>Profile9</title>
<para>



<mediaobject>
  <imageobject>
    <imagedata fileref="images/plan9mod.png" format="PNG"/>
  </imageobject>
  <imageobject>
    <imagedata fileref="images/plan9mod.eps" format="EPS"/>
  </imageobject>
</mediaobject>


The profile9 module is equal to the standard profile-HMM architecture.
The input parameter specifies the length of the module, i.e. the number of match states.
Incoming transitions connect to the an initial silent state. From this state there are transitions to an pre-model insert
state and to the first match and delete states of the model. For the purpose of local alignment with respect to
the model there are also transitions to the following match states. These may be set to zero when global alignment is
prefered. From the last match and delete states there are transitions to a silent state at the end. There are also transintions
from the previous match states to this state for the purpose of local alignments. These are set to zero when global alignment
is used. This silent end-state has a transition to an insert state for the latter unaligned part of a sequence and to the next
module.
Transition parameters are set automatically to default values at this stage. For a profile HMM, they may be updated
using the <application>opt_prfhmm</application>  program.
</para>
</sect5>


  <sect5 id="profile_7_module">
    <title>Profile7 module</title>
<para>


<mediaobject>
  <imageobject>
    <imagedata fileref="images/plan7mod.png" format="PNG"/>
  </imageobject>
  <imageobject>
    <imagedata fileref="images/plan7mod.eps" format="EPS"/>
  </imageobject>
</mediaobject>


The profile7 module is equal to the profile9 module in every way, except for it not having any delete  &#x027f9; insert or
 insert  &#x027f9; delete transitions.

</para>
</sect5>


  <sect5 id="u_turn_module">
    <title>U-turn module</title>
<para>


<mediaobject>
  <imageobject>
    <imagedata fileref="images/uturnmod.png" format="PNG"/>
  </imageobject>
  <imageobject>
    <imagedata fileref="images/uturnmod.eps" format="EPS"/>
  </imageobject>
</mediaobject>


The U-turn module is a set of states that models a symbol sequence of arbitrary length. The "bottom
of the U" is represented by a single node with a transition to itself, while the states in the leg leading to
the bottom of the U have transitions both to the next state towards the bottom and to two states on the other leg of the U.
The states in the leg leading out of the U have a transition to the next state on the path out from the U.
At initialization all transition probabilities are set to 0.5. The <varname>length</varname> parameter defines the number of
states in each leg.

</para>
</sect5>


  <sect5 id="highway_module">
    <title>Highway module</title>
<para>


<mediaobject>
  <imageobject>
    <imagedata fileref="images/highwaymod.png" format="PNG"/>
  </imageobject>
  <imageobject>
    <imagedata fileref="images/highwaymod.eps" format="EPS"/>
  </imageobject>
</mediaobject>


A highway module models a symbol sequence of arbitrary but fixed length. The <varname>length</varname> parameter defines this length.



</para>
</sect5>

</sect4>















  <sect4 id="basic_modhmmc">
    <title>basic modhmmc</title>
<para>



<variablelist 
><varlistentry><term 
>
Name of HMM? </term><listitem>
     <!--l. 9--><para>Specify file name of HMM. <application>modhmmc</application> will add '.hmg'
     </para></listitem></varlistentry><varlistentry><term 
>
Nr of alphabets (1-4)? </term><listitem>
     <!--l. 10--><para>Specify number of alphabets
     </para></listitem></varlistentry><varlistentry><term 
>
Alphabet 1: </term><listitem>
     <!--l. 11--><para>The alphabet of an HMM is specified as a set of letters, where each
     letter is a word of up to 4 characters. etc
     </para></listitem></varlistentry><varlistentry><term 
>
Alphabet 2: </term><listitem>
     <!--l. 13--><para>If multiple alphabets exist, these must also be specified
     </para></listitem></varlistentry><varlistentry><term 
>
Alphabet 3: </term><listitem>
     <!--l. 14--><para>...
     </para></listitem></varlistentry><varlistentry><term 
>
Alphabet 4: </term><listitem>
     <!--l. 15--><para>...

     </para></listitem></varlistentry><varlistentry><term 
>
Start node: </term><listitem>
     <!--l. 16--><para>Specify name of start state
     </para></listitem></varlistentry><varlistentry><term 
>
Module x: </term><listitem>
     <!--l. 17--><para>Specify module type, name and label
     </para></listitem></varlistentry><varlistentry><term 
>
End node: </term><listitem>
     <!--l. 18--><para>Specify name of end state
     </para></listitem></varlistentry><varlistentry><term 
>
Specify interconnectivity </term><listitem>
     <!--l. 19--><para>
     </para></listitem></varlistentry><varlistentry><term 
>
Connection from a to: </term><listitem>
     <!--l. 20--><para>Specify transitions from state a
     </para></listitem></varlistentry><varlistentry><term 
>
Specify emission distribution groups: </term><listitem>
     <!--l. 21--><para>Tie the emission probabilities of (all states of) the specified modules to
     each other. This means that when updating the value of the emission
     probabilities during training, all states in a distribution group will be
     treated as one state.
     </para></listitem></varlistentry><varlistentry><term 
>
Specify transition distribution groups: </term><listitem>
     <!--l. 24--><para>Tie the transition probabilities of the specified modules to each other.
     Modules must be identical (same type and size) for this to work. Tying
     two modules means that the transitions in the two modules which
     corresponds to each other will be updated as the one transition when
     transition probabilities are updated during training.
     </para></listitem></varlistentry><varlistentry><term 
>
Specify initial emission probabilities ... </term><listitem>
     <!--l. 27--><para>Initialize the emission probabilities and tie the states to specific prior
     distribution files (see <xref linkend="distribution_files"/>
). Also specify the weight of the prior against
     the observations, default is 1.0. This is done for each alphabet.
     </para></listitem></varlistentry><varlistentry><term 
>
Specify initial transition probabilities ... </term><listitem>
     <!--l. 31--><para>Same  as  above,  but  for  transitions.  (prior  files  and  weights  not
     implemented in training and scoring algorithms yet).</para></listitem></varlistentry></variablelist>






</para>


  </sect4>

  <sect4 id="modhmmc_examples">
    <title>Examples</title>
<para>

No examples yet
<!--
<example id="modhmmc_asdfasdf"><title>fastdist withasdfsdf</title><para>
We use the file described in <xref linkend="seq.phylip_multialignment"/> as input file. 
The file has two datasets so we pass the option <userinput>-r 2</userinput> to <application>fastdist</application>. Per default the output is given in XML format 

<programlisting><![CDATA[
[user@saturn ~]$ ]]><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/commands/fastdist_seq.phylip_multialignment">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
<xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_BINARY_DIR}/xincluded_files/fastdist_seq.phylip_multialignment">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
</programlisting>


</para>
</example>

-->


      </para>

    </sect4>

    </sect3>


  <sect3 id="modhmmt">
    <title>modhmmt</title>
<para>
The <application>modhmmt</application> program is is used for parameter optimization. It implements the regular Baum-Welch training algorithm as well as
Conditional maximum likelihood (CML) training, both for single sequences, multiple sequence alignments and sequence profiles.
</para>



  <sect4 id="modhmmt_command_line_options">
    <title>Command line options</title>


<para>

Type <userinput>modhmmt --help</userinput> to see the command line options

<programlisting><![CDATA[
[user@saturn ~]$ ]]><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/commands/modhmmt_help">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
<xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_BINARY_DIR}/xincluded_files/modhmmt_help">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
</programlisting>

</para></sect4>



  <sect4 id="modhmmt_examples">
    <title>Examples</title>
<para>

no examples yet
<!--
<example id="modhmmt_phylip.dm"><title>modhmmt with input file in Phylip distance matrix format</title><para>
We use the file described in <xref linkend="dm.phylip_dm"/> as input file. The file has two datasets so we pass the option <userinput>-r 2</userinput> to <application>modhmmt</application>. Per default the output is given in the "modhmm count tree XML format" ( -O xml ).

<programlisting><![CDATA[
[user@saturn ~]$ ]]><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/commands/modhmmt_dm.phylip_dm">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
<xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_BINARY_DIR}/xincluded_files/modhmmt_dm.phylip_dm">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
</programlisting>
</para>
</example>

-->
</para>
</sect4>
    </sect3>



  <sect3 id="modhmms">
    <title>modhmms</title>
<para>
The <application>modhmms</application> program is used for scoring sequences, multiple sequence alignments and sequence profiles against an hmm. The algorithms implemented
includ forward, Viterbi and 1-best. Output can be either a log-likelihood/logodds/reverse score, the (approximatley) most probable labeling of a sequence
or the most probale state path.
</para>



  <sect4 id="modhmms_command_line_options">
    <title>Command line options</title>


<para>

Type <userinput>modhmms --help</userinput> to see the command line options

<programlisting><![CDATA[
[user@saturn ~]$ ]]><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/commands/modhmms_help">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
<xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_BINARY_DIR}/xincluded_files/modhmms_help">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
</programlisting>

</para></sect4>


  <sect4 id="modhmms_examples">
    <title>Examples</title>
<para>

no examples yet
<!--
<example id="modhmms_phylip.dm"><title>modhmms with input file in Phylip distance matrix format</title><para>
We use the file described in <xref linkend="dm.phylip_dm"/> as input file. The file has two datasets so we pass the option <userinput>-r 2</userinput> to <application>modhmms</application>. Per default the output is given in the "modhmm count tree XML format" ( -O xml ).

<programlisting><![CDATA[
[user@saturn ~]$ ]]><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/commands/modhmms_dm.phylip_dm">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
<xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_BINARY_DIR}/xincluded_files/modhmms_dm.phylip_dm">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
</programlisting>
</para>
</example>

-->
</para>
</sect4>

    </sect3>


  <sect3 id="modseqalign">
    <title>modseqalign</title>
<para>

The <application>modseqalign</application> program aligns two sequences/multiple alignments/profiles using their most proable state paths through a given hmm.
</para>



  <sect4 id="modseqalign_command_line_options">
    <title>Command line options</title>


<para>

Type <userinput>modseqalign --help</userinput> to see the command line options

<programlisting><![CDATA[
[user@saturn ~]$ ]]><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/commands/modseqalign_help">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
<xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_BINARY_DIR}/xincluded_files/modseqalign_help">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
</programlisting>

</para></sect4>




  <sect4 id="modseqalign_examples">
    <title>Examples</title>
<para>

no examples yet
<!--
<example id="modseqalign_phylip.dm"><title>modseqalign with input file in Phylip distance matrix format</title><para>
We use the file described in <xref linkend="dm.phylip_dm"/> as input file. The file has two datasets so we pass the option <userinput>-r 2</userinput> to <application>modseqalign</application>. Per default the output is given in the "modhmm count tree XML format" ( -O xml ).

<programlisting><![CDATA[
[user@saturn ~]$ ]]><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/commands/modseqalign_dm.phylip_dm">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
<xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_BINARY_DIR}/xincluded_files/modseqalign_dm.phylip_dm">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
</programlisting>
</para>
</example>

-->
</para>
</sect4>

    </sect3>


  <sect3 id="add_alphabet">
    <title>add_alphabet</title>
<para>
The <application>add_alphabet</application> program is used for scoring sequences, multiple sequence alignments and sequence profiles against an hmm. The algorithms implemented
includ forward, Viterbi and 1-best. Output can be either a log-likelihood/logodds/reverse score, the (approximatley) most probable labeling of a sequence
or the most probale state path.
</para>



  <sect4 id="add_alphabet_command_line_options">
    <title>Command line options</title>


<para>

Type <userinput>add_alphabet --help</userinput> to see the command line options

<programlisting><![CDATA[
[user@saturn ~]$ ]]><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/commands/add_alphabet_help">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
<xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_BINARY_DIR}/xincluded_files/add_alphabet_help">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
</programlisting>

</para></sect4>


  <sect4 id="add_alphabet_examples">
    <title>Examples</title>
<para>

no examples yet
<!--
<example id="add_alphabet_phylip.dm"><title>add_alphabet with input file in Phylip distance matrix format</title><para>
We use the file described in <xref linkend="dm.phylip_dm"/> as input file. The file has two datasets so we pass the option <userinput>-r 2</userinput> to <application>add_alphabet</application>. Per default the output is given in the "modhmm count tree XML format" ( -O xml ).

<programlisting><![CDATA[
[user@saturn ~]$ ]]><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/commands/add_alphabet_dm.phylip_dm">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
<xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_BINARY_DIR}/xincluded_files/add_alphabet_dm.phylip_dm">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include>
</programlisting>
</para>
</example>

-->
</para>
</sect4>

    </sect3>



    </sect2>

  <sect2 id="file_formats">
    <title>File formats</title>
    <para>      





There are a few different types of input files accociated with the modhmm package, mainly files
for sequences and various utility files. There are 4 different sequence file formats (fasta, single, multi and profile).
The profile format is compulsory to contain labels, while for the other formats labels are optional. Each sequence file
may contain only one sequence, but a sequence may consist of up to 4 different parallel alphabets.




</para>

  <sect3 id="file_formats_for_discrete_alphabets">
    <title>File formats for discrete alphabets</title>
      



  <sect4 id="modhmm_fa_format">
    <title>Modhmm fa format</title>
<para>
The standard fasta sequence file format. Fasta sequences can only consist of one alphabet. 


<example id="seq.modhmm-fa"><title>seq.modhmm-fa, an example file in the Modhmm fa format</title><para>
The example file <ulink url="example_files/seq.modhmm-fa">seq.modhmm-fa</ulink> contains


<programlisting><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/example_files/seq.modhmm-fa">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include></programlisting>

</para></example>

Labels can also be used in this format  ( see <xref linkend="labels"/> ).

</para></sect4>



  <sect4 id="modhmm_s_discrete_format">
    <title>Modhmm s discrete format</title>
<para>


The <emphasis>Modhmm s discrete format</emphasis> is a single sequence file format specific for modhmm, designed to allow for an alphabet which allows letters with more than one character.
'&lt;' is the signal for starting a sequence, '&gt;' is the signal for ending it and ';' is the end-of-letter marker. Note that
this marker must also be present after the last letter in a sequences, i.e. right before the '&gt;' sign. A sequence
may occupy more than one line and can be of arbitrary length. Anything written outside of a &lt;&gt; marking will be ignored.

<example id="seq.modhmm-s-discrete"><title>seq.modhmm-s-discrete, an example file in the Modhmm s discrete format</title><para>
The example file <ulink url="example_files/seq.modhmm-s-discrete">seq.modhmm-s-discrete</ulink> contains
<programlisting><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/example_files/seq.modhmm-s-discrete">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include></programlisting>
</para></example>


The multiple alphabet version of this format has the same basic structure as the single alphabet version. Each sequence
consists of 2-4 subsequences which begin and end with '&lt;' and '&gt;' respectively. The number of '&lt;' characters in the
beginning is used to indicate which alphabet a subsequence
belongs to. The subsequences must be in numerical order, and they must be of equal lengths.





<example id="seq-multiple-alphabet.modhmm-s-discrete"><title>seq-multiple-alphabet.modhmm-s-discrete, a multiple-alphabet example file in the Modhmm s discrete format</title><para>
The example file <ulink url="example_files/seq-multiple-alphabet.modhmm-s-discrete">seq-multiple-alphabet.modhmm-s-discrete</ulink> contains
<programlisting><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/example_files/seq-multiple-alphabet.modhmm-s-discrete">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include></programlisting>
</para></example>


Labels can also be used in this format  ( see <xref linkend="labels"/> ).
</para></sect4>




  <sect4 id="modhmm_msa_discrete_format">
    <title>Modhmm msa discrete format</title>
<para>
The <emphasis>Modhmm msa discrete format</emphasis> is the same format as the <emphasis>Modhmm s discrete format</emphasis> ( <xref linkend="modhmm_s_discrete_format"/> ) but for multiple sequence alignments. This means that one msa sequence file
contains a set of sequences aligned to each other. For this purpose, each sequence must be situated in one line only.
'-', '_', ' ' and '.' can all be used to signify a gap. All sequences in a msa sequence file are regarded as belonging to
the same multiple sequence alignment. Each line must be as long as or shorter than the line containing the template sequence
when such a sequence is used. The starting position of an alignment is calculated from the first letter after '&lt;',
which means that it is not possible to skip letters in the begining of a sequence. These must be represented as gaps.
However, it is possible to do this at the end of a sequence.





<example id="seq.modhmm-msa-discrete"><title>seq.modhmm-msa-discrete, an example file in Modhmm msa discrete format</title><para>
The example file <ulink url="example_files/seq.modhmm-msa-discrete">seq.modhmm-msa-discrete</ulink> contains



<programlisting><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/example_files/seq.modhmm-msa-discrete">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include></programlisting>

It uses a single alphabet.
</para>
</example>




This however will not work:

<programlisting><![CDATA[<A;B;C;D;E;F;G;H;>
<A;-;C;D;E;P;G;H;>
<A;B;-;-;E; ; ; ;>
    <C;D;-;F>
]]></programlisting>
'C' in the fourth line will be aligned with the 'A'-column, 'D' with the 'B'-column and so on.

<example id="seq-multiple-alphabet.modhmm-msa-discrete"><title>seq-multiple-alphabet.modhmm-msa-discrete, an example file in the Modhmm msa discrete format</title><para>
The example file <ulink url="example_files/seq-multiple-alphabet.modhmm-msa-discrete">seq-multiple-alphabet.modhmm-msa-discrete</ulink> contains

<programlisting><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/example_files/seq-multiple-alphabet.modhmm-msa-discrete">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include></programlisting>

It uses multiple alphabets.
</para>
</example>

Labels can also be used in this format  ( see <xref linkend="labels"/> ).
</para></sect4>




  <sect4 id="labels">
    <title>Labels</title>
<para>

A line containing a labeling of the sequence may be added to sequences in fa, s or msa format. The
labeling is represented as a line at the end of the file (enclosed in by '/' characters) with one label character for each sequence position.
For a fasta sequence:

<programlisting><![CDATA[>seqname1
ACDEFGHIKLMNPQRSTVWY
/..................../
]]></programlisting>

For an s sequence:

<programlisting><![CDATA[>seqname1
<A;C;D;E;F;G;H;I;K;L;M;N;P;Q;R;S;T;V;W;Y;>
/..................../
]]></programlisting>

For an msa sequence:

<programlisting><![CDATA[>seqname1
<A;C;D;E;F;G;H;I;K;L;M;N;P;Q;R;S;T;V;W;Y;>
<A;C;D;T;F;T;H;I;K;-;-;N;P;Q;R;S;T;-;W;Y;>
<-;-;-;-;-;-;-;-;-;-;-;N;P;Q;R;S;T;-;Y;Y;>
/..................../
]]></programlisting>


A label must only consist of one character. The character '.' is predefined as representing any label, which means that a
sequence letter labeled with '.' will match any state label.

For multiple alphabets the principle is the same as for unlabeled alignment sequences, and the labels are placed at the bottom
(only once, since the labeling is the same over all alphabets).




</para></sect4>
</sect3>


  <sect3 id="file_formats_for_continuous_alphabets">
    <title>File formats for continuous alphabets</title>







  <sect4 id="modhmm_s_continuous_format">
    <title>Modhmm s continuous format</title>
<para>


The <emphasis>Modhmm s continuous format</emphasis>  is the continuos version of the <emphasis>Modhmm s discrete format</emphasis>. All lines start and end with '#' and '+' respectively
to mark that the alphabet is continuous. Otherwise, there is no difference from the discrete std format.



<example id="seq.modhmm-s-continuous"><title>seq.modhmm-s-continuous, an example file in the Modhmm s continuous format</title><para>
The example file <ulink url="example_files/seq.modhmm-s-continuous">seq.modhmm-s-continuous</ulink> contains
<programlisting><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/example_files/seq.modhmm-s-continuous">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include></programlisting>
</para></example>

The multiple alphabet version of this format has the same basic structure as the single alphabet version. Each sequence
consists of 2-4 subsequences which begin and end with '#' and '+' respectively. The number of '#' characters in the
beginning is used to indicate which alphabet a subsequence
belongs to. The subsequences must be in numerical order, and they must be of equal lengths. It is also possible to
mix discrete and continuous alphabets. This is simply done by using the indicators '#'- '+' and '&lt;' and '&gt;' respectively.





<example id="seq-multiple-alphabet.modhmm-s-continuous"><title>seq-multiple-alphabet.modhmm-s-continuous, an example file in the Modhmm s continuous format</title><para>
The example file <ulink url="example_files/seq-multiple-alphabet.modhmm-s-continuous">seq-multiple-alphabet.modhmm-s-continuous</ulink> contains
<programlisting><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/example_files/seq-multiple-alphabet.modhmm-s-continuous">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include></programlisting>
</para></example>


Labeling sequences from a continuous alphabet works in exactly the same way as for sequences from a discrete alphabet.


</para>

</sect4>




  <sect4 id="modhmm_msa_continuous_format">
    <title>Modhmm msa continuous format</title>
<para>


Alignments of sequences from continuous alphabets do not exist. However,
when using multiple alphabets, single sequences from continuous alphabets may be used together
with alignments of ones from discrete alphabets. Here, 'X' is a wildcard letter which stands for anything
(emission probability 1.0).



<example id="seq.modhmm-msa-continuous"><title>seq.modhmm-msa-continuous, an example file in Modhmm msa continuous format</title><para>
The example file <ulink url="example_files/seq.modhmm-msa-continuous">seq.modhmm-msa-continuous</ulink> contains



<programlisting><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/example_files/seq.modhmm-msa-continuous">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include></programlisting>


</para>
</example>
</para></sect4>


</sect3>

  <sect3 id="modhmm_prf_format">
    <title>Modhmm prf format</title>
<para>

The <emphasis>Modhmm prf format</emphasis> is slightly different form the other formats as it does not contain letters, but rather frequencies at each position. Labels are
compulsory for this sequence format.

<example id="seq.modhmm-prf"><title>seq.modhmm-prf, an example file in Modhmm msa prf format</title><para>
The example file <ulink url="example_files/seq.modhmm-prf">seq.modhmm-prf</ulink> contains



<programlisting><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/example_files/seq.modhmm-prf">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include></programlisting>


</para>
</example>
</para></sect3>





  <sect3 id="modhmm_replacement_letter_format">
    <title>Modhmm replacement letter format</title>
<para>

The <emphasis>Modhmm replacement letter format</emphasis><![CDATA[ is fairly strict. All lines beginning with '#' (as the very first character), and all
empty lines are disregarded. Apart from this, the first line must be the number of alphabets, the second line must be a
number giving the size of the first alphabet.
The third line must be the characters of the alphabet, separated by and ending with ';'. The
fourth line is the number of the replacement letters.
Next, depending on the number in the fourth line, there are a number of lines with the specification of the replacement letters for that alphabet.
Then the next line contains the number giving the size of the second alphabet if it exists, etc.]]></para><para><![CDATA[

The specification lines have the following layout: <alphabet letter> = <substitution> <substitution>
<substitution> ... . The substitutions have the form: <replacement letter>:<probabilityshare>.
]]>

<example id="file.modhmm-replacement-letter"><title>file.modhmm-replacement-letter, an example file in Modhmm replacement letter format</title><para>
The example file <ulink url="example_files/file.modhmm-replacement-letter">file.modhmm-replacement-letter</ulink> contains



<programlisting><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/example_files/file.modhmm-replacement-letter">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include></programlisting>


</para>
</example>
</para></sect3>




  <sect3 id="modhmm_subst_mtx_format">
    <title>Modhmm subst mtx format</title>
<para>

The <emphasis>Modhmm subst mtx format</emphasis> specifies a substitution matrix as defined in <xref linkend="emission_scoring_methods"/><![CDATA[.
All lines beginning with '#' (as the very first character), and all
empty lines are disregarded. Apart from this, the first line must be a number giving the size of the alphabet of the associated hmm.
The second line must be the characters of the alphabet, separated by and ending with ';'.
The following lines specify the matrix itself. Each line has the following layout: <alphabet letter> = <score alphabet letter 1>
<score alphabet letter 2> ... .
There is one line for each letter of the alphabet, not in any particular order, and one line at the end that describes
how to interpret letters in the sequences not in the alphabet. This line is only used under certain circumstances in the
SMDP and SMDPP scoring methods. There is one column for each letter of the alphabet and the columns must be in the
same order as the alphabet specification of the hmm.
 ]]>
  

<example id="file.modhmm-subst-mtx"><title>file.modhmm-subst-mtx, an example file in Modhmm subst mtx format</title><para>
The example file <ulink url="example_files/file.modhmm-subst-mtx">file.modhmm-subst-mtx</ulink> contains
<programlisting><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/example_files/file.modhmm-subst-mtx">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include></programlisting>





</para>
</example>
</para></sect3>





  <sect3 id="modhmm_frequency_format">
    <title>Modhmm frequency format</title>
<para>

The <emphasis>Modhmm frequency format</emphasis> specifies the background frequencies of the alphabet letters in some chosen population
as defined in <xref linkend="emission_scoring_methods"/><![CDATA[.
All lines beginning with '#' (as the very first character), and all
empty lines are disregarded. Apart from this, the first line must be a number giving the size of the alphabet of the associated hmm.
The following lines specify the frequencies.
There is one line for each letter of the alphabet, in the same order as the alphabet is listed in the hmm specification.


]]>
<example id="file.modhmm-frequency"><title>file.modhmm-frequency, an example file in Modhmm frequency format</title><para>
The example file <ulink url="example_files/file.modhmm-frequency">file.modhmm-frequency</ulink> contains
<programlisting><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/example_files/file.modhmm-frequency">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include></programlisting>





</para>
</example>




Note that there are no multiple alphabet frequency files. These are specified separately for each alphabet.


  
</para></sect3>



  <sect3 id="modhmm_prior_format">
    <title>Modhmm prior format</title>
<para>

All lines beginning with '#' (as the very first character), and all
empty lines are disregarded. Apart from this, the first line should contain the number of mixture components.
Then, for each component there is a line containing the component's probability value
and a line containing the component values for each alphabet letter. Note that the order of the component values
must match the alphabet order in the hmm specification, i.e. the first component value of a line will be associated
with the first alphabet letter in the hmm specification, etc.


<example id="file.modhmm-prior"><title>file.modhmm-prior, an example file in Modhmm prior format</title><para>
The example file <ulink url="example_files/file.modhmm-prior">file.modhmm-prior</ulink> contains
<programlisting><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/example_files/file.modhmm-prior">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include></programlisting>





</para>
</example>


Note that there are no multiple alphabet priorfiles. These are specified separately for each alphabet.



  
</para></sect3>




  <sect3 id="modhmm_null_model_format">
    <title>Modhmm null model format</title>
<para>
All lines beginning with '#' (as the very first character), and all
empty lines are disregarded. Apart fom this, the first line should contain the alphabet size (<varname>n</varname>). The lines
2 to <varname>n</varname>+1, should contain the emission probabilities for the letters of the alphabet, ordered the same way as
in the hmm specification. The last line should contain the loop transition probability. The null model is
built as an hmm wit one singleloop module.

<example id="file.modhmm-null-model"><title>file.modhmm-null-model, an example file in Modhmm null model format</title><para>
The example file <ulink url="example_files/file.modhmm-null-model">file.modhmm-null-model</ulink> contains
<programlisting><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/example_files/file.modhmm-null-model">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include></programlisting>



</para>
</example>

 
</para></sect3>

  <sect3 id="modhmm_name_file_format">
    <title>Modhmm name file format</title>
<para>



The files called hmmnamefile and seqnamefile in the options are in the <emphasis>Modhmm name file format</emphasis>. These files are simply files with names
of hmms/sequences inluding a relative or full path. Note that this file format is very sensitive to blank lines, blanks
added at the end of lines etc. Each line should contain only the name of the hmm/sequence file.



<example id="file.modhmm-name-file"><title>file.modhmm-name-file, an example file in Modhmm name file format</title><para>
The example file <ulink url="example_files/file.modhmm-name-file">file.modhmm-name-file</ulink> contains
<programlisting><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/example_files/file.modhmm-name-file">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include></programlisting>



</para>
</example>

 
</para></sect3>


  <sect3 id="modhmm_alphabet_format">
    <title>Modhmm alphabet format</title>
<para>
Files in the <emphasis>Modhmm alphabet format</emphasis> are used to specify the alphabet and the emission probabilities when adding an alphabet to
an already existing HMM. The first row contains the text 'ALPHABET: ' and a string containing the alphabet where the
letters are separated by a ';'. The second row conatains the text 'ALPHABET LENGTH: ' and a positive integer giving
the number of letters in the alphabet. The remaining rows contain the text 'VERTEX X: ' and a sequence of floating point
numbers (separated by a 'SPACE' character). The sequence of floating point numbers give the emission probabilities for the
alphabet letters in state X of the model. This means that there should be as many floats in the sequence of each line as there
are letters in the alphabet. The order of the probabilities should follow the order of the letters in the alphabet specification
from the first row. There should be one row for each state in the HMM, which the alphabet is to be added to. The sum of the
numbers from one row should be either 1.0 (for emitter states) or 0.0 (for silent states) for discrete alphabets.



<example id="file.modhmm-alphabet"><title>file.modhmm-alphabet, an example file in Modhmm alphabet format</title><para>
The example file <ulink url="example_files/file.modhmm-alphabet">file.modhmm-alphabet</ulink> contains
<programlisting><xi:include xmlns:xi="http://www.w3.org/2001/XInclude" parse="text" encoding="UTF-8" href="${CMAKE_CURRENT_SOURCE_DIR}/example_files/file.modhmm-alphabet">
<xi:fallback>
   couldn't xinclude file
</xi:fallback>
</xi:include></programlisting>



</para>
</example>

 
</para></sect3>












  </sect2>


  </sect1>


  <sect1 id="functionality">
    <title>Functionality</title>

  <sect2 id="emission_scoring_methods">
    <title>Emission scoring methods</title>
    <para>


There are 7 different ways to calculate the score (probability) for a multiple sequence alignment column to be emitted
in a certain state. Let us call this score  
<inlineequation>  
<mml:math>
<mml:mrow>
 <mml:msup>
   <mml:mi>e</mml:mi>
   <mml:mi>msa</mml:mi>
 </mml:msup>
</mml:mrow>
  </mml:math>
 </inlineequation>
  .

</para>

  <sect3 id="emission_scoring_method_dp">
    <title>DP</title>
    <para>







The first scoring method, DP,  is the dot product


<equation>
  <mml:math>
    <mml:mstyle displaystyle="true">
    <mml:mrow>
      <mml:msubsup>
        <mml:mi>e</mml:mi>
        <mml:mi>q</mml:mi>
        <mml:mi>msa</mml:mi>
      </mml:msubsup>
      <mml:mo  stretchy="false">(</mml:mo>
      <mml:msubsup>
        <mml:mi>x</mml:mi>
        <mml:mi>k</mml:mi>
        <mml:mi>msa</mml:mi>
      </mml:msubsup>
      <mml:mo  stretchy="false">)</mml:mo>
      <mml:mo>=</mml:mo>
      <mml:munder>
        <mml:mrow>
          <mml:mo>
&sum;    </mml:mo>
        </mml:mrow>
        <mml:mi>j</mml:mi>
      </mml:munder>
      <mml:msub>
        <mml:mi>e</mml:mi>
        <mml:mi>q</mml:mi>
      </mml:msub>
      <mml:mo  stretchy="false">(</mml:mo>
      <mml:msub>
        <mml:mi>&sigma;</mml:mi>
        <mml:mi>j</mml:mi>
      </mml:msub>
      <mml:mo  stretchy="false">)</mml:mo>
      <mml:mo>*</mml:mo>
      <mml:msubsup>
        <mml:mi>x</mml:mi>
        <mml:mi>k</mml:mi>
        <mml:mi>msa</mml:mi>
      </mml:msubsup>
      <mml:mo  stretchy="false">(</mml:mo>
      <mml:msub>
        <mml:mi>&sigma;</mml:mi>
        <mml:mi>j</mml:mi>
      </mml:msub>
      <mml:mo  stretchy="false">)</mml:mo>
    </mml:mrow>
    </mml:mstyle>
  </mml:math>
</equation>




Here,
<inlineequation>  
 <mml:math>
  <mml:msub>
   <mml:mi>&sigma;</mml:mi>
   <mml:mi>j</mml:mi>
  </mml:msub>
 </mml:math>
</inlineequation>

 represents alphabet symbol
<inlineequation>  
<mml:math>
   <mml:mi>j</mml:mi>
  </mml:math>
 </inlineequation>, 
<inlineequation>  
<mml:math>
   <mml:mi>q</mml:mi>
  </mml:math>
 </inlineequation>

 is s state and 

<inlineequation>  
<mml:math>
 <mml:msubsup>
   <mml:mi>x</mml:mi>

   <mml:mi>k</mml:mi>
   <mml:mi>msa</mml:mi>
 </mml:msubsup>

  </mml:math>
 </inlineequation>

 is sequence column <inlineequation>  
<mml:math>

   <mml:mi>k</mml:mi>


  </mml:math>
 </inlineequation>
.
This formula is easy to use and only increases the time complexity of the forward and backward algorithms by a factor
the size of the alphabet. The probabilistic interpretation of this score is that it represents the probability that the same
letter is drawn if drawing independently from both the state- and the profile- distributions.
</para>
</sect3>




  <sect3 id="emission_scoring_method_gm">
    <title>GM</title>
<para>
The second implementation of <inlineequation>  
<mml:math>
<mml:mrow>
 <mml:msup>
   <mml:mi>e</mml:mi>
   <mml:mi>msa</mml:mi>
 </mml:msup>
</mml:mrow>
  </mml:math>
 </inlineequation>
, GM, calculates the total probability of scoring all sequences of the profile
as single sequences but forced to take the same path through the model, normalized on the number of sequences in
the alignment






<equation>
  <mml:math>

    <mml:mstyle displaystyle="true">
    <mml:mrow>
      <mml:msubsup>
        <mml:mi>e</mml:mi>
        <mml:mi>q</mml:mi>
        <mml:mi>msa</mml:mi>
      </mml:msubsup>
      <mml:mo  stretchy="false">(</mml:mo>
      <mml:msubsup>
        <mml:mi>x</mml:mi>
        <mml:mi>k</mml:mi>
        <mml:mi>msa</mml:mi>
      </mml:msubsup>
      <mml:mo  stretchy="false">)</mml:mo>
      <mml:mo>=</mml:mo>
      <mml:munder>
        <mml:mrow>
          <mml:mo>
&prod;    </mml:mo>
        </mml:mrow>
        <mml:mi>j</mml:mi>
      </mml:munder>
      <mml:msup>
        <mml:mrow>
          <mml:msub>
            <mml:mi>e</mml:mi>
            <mml:mi>q</mml:mi>
          </mml:msub>
          <mml:mo  stretchy="false">(</mml:mo>
          <mml:msub>
            <mml:mi>&sigma;</mml:mi>
            <mml:mi>j</mml:mi>
          </mml:msub>
          <mml:mo  stretchy="false">)</mml:mo>
        </mml:mrow>
        <mml:mrow>
          <mml:msubsup>
            <mml:mi>x</mml:mi>
            <mml:mi>k</mml:mi>
            <mml:mi>msa</mml:mi>
          </mml:msubsup>
          <mml:mo  stretchy="false">(</mml:mo>
          <mml:msub>
            <mml:mi>&sigma;</mml:mi>
            <mml:mi>j</mml:mi>
          </mml:msub>
          <mml:mo  stretchy="false">)</mml:mo>
        </mml:mrow>
      </mml:msup>
    </mml:mrow>
    </mml:mstyle>
  </mml:math>
</equation>


As for the scoring method above, this formula only increases the time complexity of the forward and backward algorithms
with a factor the size of the alphabet.


</para>

</sect3>
  <sect3 id="emission_scoring_method_picasso">
    <title>PICASSO</title>

<para>
The third scoring method, PICASSO, is a scoring method originally developed for non-HMM-based profile-profile
comparison. The probabilistic interpretaion of this score, is that it (as GM) calculates the geometric mean, but
an emission probability is not a pure frequency, but a frequency which is weighted according to how common this
letter is in a given background distribution. In the symmetric version, the profile columns are weighted in a
similar fashion. The respective columns are then normalized to sum to 1. In PICASSO, the emission score is calculated as


<equation>
  <mml:math>
    <mml:mstyle displaystyle="true">
    <mml:mrow>
      <mml:msubsup>
        <mml:mi>e</mml:mi>
        <mml:mi>q</mml:mi>
        <mml:mi>msa</mml:mi>
      </mml:msubsup>
      <mml:mo  stretchy="false">(</mml:mo>
      <mml:msubsup>
        <mml:mi>x</mml:mi>
        <mml:mi>k</mml:mi>
        <mml:mi>msa</mml:mi>
      </mml:msubsup>
      <mml:mo  stretchy="false">)</mml:mo>
      <mml:mo>=</mml:mo>
      <mml:munderover>
          <mml:mo>
&prod;    </mml:mo>
        <mml:mrow>
         <mml:mi>j</mml:mi>
          <mml:mo>=</mml:mo>
          <mml:mn>1</mml:mn>
        </mml:mrow>
        <mml:mi>Z</mml:mi>
      </mml:munderover>
      <mml:msup>

        <mml:mrow>
       <mml:mo     >(</mml:mo>
        <mml:mfrac>
        <mml:mrow>
          <mml:msub>
            <mml:mi>e</mml:mi>
            <mml:mi>q</mml:mi>
          </mml:msub>
          <mml:mo stretchy="false">(</mml:mo>
          <mml:msub>
            <mml:mi>&sigma;</mml:mi>
            <mml:mi>j</mml:mi>
          </mml:msub>
          <mml:mo stretchy="false">)</mml:mo>
        </mml:mrow>
         <mml:msub>
            <mml:mi>freq</mml:mi>
          <mml:msub>
            <mml:mi>&sigma;</mml:mi>
            <mml:mi>j</mml:mi>
          </mml:msub>
          </mml:msub>
        </mml:mfrac>
       <mml:mo  >)</mml:mo>

        </mml:mrow>
        <mml:mrow>
          <mml:msubsup>
            <mml:mi>x</mml:mi>
            <mml:mi>k</mml:mi>
            <mml:mi>msa</mml:mi>
          </mml:msubsup>
          <mml:mo  stretchy="false">(</mml:mo>
          <mml:msub>
            <mml:mi>&sigma;</mml:mi>
            <mml:mi>j</mml:mi>
          </mml:msub>
          <mml:mo  stretchy="false" >)</mml:mo>
        </mml:mrow>
      </mml:msup>
    </mml:mrow>
    </mml:mstyle>
  </mml:math>
</equation>





for the original version 
<biblioref linkend="heger:2003"/>
 or as 


<equation>
  <mml:math>
    <mml:mstyle displaystyle="true">
    <mml:mrow>
      <mml:msubsup>
        <mml:mi>e</mml:mi>
        <mml:mi>q</mml:mi>
        <mml:mi>msa</mml:mi>
      </mml:msubsup>
      <mml:mo  stretchy="false">(</mml:mo>
      <mml:msubsup>
        <mml:mi>x</mml:mi>
        <mml:mi>k</mml:mi>
        <mml:mi>msa</mml:mi>
      </mml:msubsup>
      <mml:mo  stretchy="false">)</mml:mo>
      <mml:mo>=</mml:mo>
      <mml:munderover>
          <mml:mo>
&prod;    </mml:mo>
        <mml:mrow>
         <mml:mi>j</mml:mi>
          <mml:mo>=</mml:mo>
          <mml:mn>1</mml:mn>
        </mml:mrow>
        <mml:mi>Z</mml:mi>
      </mml:munderover>
      <mml:msup>

        <mml:mrow>
       <mml:mo     >(</mml:mo>
        <mml:mfrac>
        <mml:mrow>
          <mml:msub>
            <mml:mi>e</mml:mi>
            <mml:mi>q</mml:mi>
          </mml:msub>
          <mml:mo stretchy="false">(</mml:mo>
          <mml:msub>
            <mml:mi>&sigma;</mml:mi>
            <mml:mi>j</mml:mi>
          </mml:msub>
          <mml:mo stretchy="false">)</mml:mo>
        </mml:mrow>
         <mml:msub>
            <mml:mi>freq</mml:mi>
          <mml:msub>
            <mml:mi>&sigma;</mml:mi>
            <mml:mi>j</mml:mi>
          </mml:msub>
          </mml:msub>
        </mml:mfrac>
       <mml:mo  >)</mml:mo>

        </mml:mrow>
        <mml:mrow>
          <mml:msubsup>
            <mml:mi>x</mml:mi>
            <mml:mi>k</mml:mi>
            <mml:mi>msa</mml:mi>
          </mml:msubsup>
          <mml:mo  stretchy="false">(</mml:mo>
          <mml:msub>
            <mml:mi>&sigma;</mml:mi>
            <mml:mi>j</mml:mi>
          </mml:msub>
          <mml:mo  stretchy="false" >)</mml:mo>
        </mml:mrow>
      </mml:msup>

      <mml:mo   >*</mml:mo>
      <mml:msubsup>
        <mml:mi>e</mml:mi>
        <mml:mi>q</mml:mi>
        <mml:mi>msa</mml:mi>
      </mml:msubsup>
      <mml:mo  stretchy="false">(</mml:mo>
      <mml:msubsup>
        <mml:mi>x</mml:mi>
        <mml:mi>k</mml:mi>
        <mml:mi>msa</mml:mi>
      </mml:msubsup>
      <mml:mo  stretchy="false">)</mml:mo>

<mml:mo   >=</mml:mo>


      <mml:munderover>
          <mml:mo>
&prod;    </mml:mo>
        <mml:mrow>
         <mml:mi>j</mml:mi>
          <mml:mo>=</mml:mo>
          <mml:mn>1</mml:mn>
        </mml:mrow>
        <mml:mi>Z</mml:mi>
      </mml:munderover>
      <mml:msup>

        <mml:mrow>
       <mml:mo     >(</mml:mo>
        <mml:mfrac>

        <mml:mrow>
     <mml:msubsup>
        <mml:mi>x</mml:mi>
        <mml:mi>k</mml:mi>
        <mml:mi>msa</mml:mi>
      </mml:msubsup>
      <mml:mo  stretchy="false">(</mml:mo>
      <mml:msub>
        <mml:mi>&sigma;</mml:mi>
        <mml:mi>j</mml:mi>
      </mml:msub>
      <mml:mo  stretchy="false">)</mml:mo>
        </mml:mrow>
         <mml:msub>
            <mml:mi>freq</mml:mi>
          <mml:msub>
            <mml:mi>&sigma;</mml:mi>
            <mml:mi>j</mml:mi>
          </mml:msub>
          </mml:msub>
        </mml:mfrac>
       <mml:mo  >)</mml:mo>

        </mml:mrow>
        <mml:mrow>
         <mml:msub>
            <mml:mi>e</mml:mi>
            <mml:mi>q</mml:mi>
          </mml:msub>
          <mml:mo stretchy="false">(</mml:mo>
          <mml:msub>
            <mml:mi>&sigma;</mml:mi>
            <mml:mi>j</mml:mi>
          </mml:msub>
          <mml:mo stretchy="false">)</mml:mo>
        </mml:mrow>
      </mml:msup>

    </mml:mrow>
    </mml:mstyle>
  </mml:math>
</equation>



for the symmetric version 
<biblioref linkend="mittelman:2003"/>.


Here,
<inlineequation>  
<mml:math>
         <mml:msub>
            <mml:mi>freq</mml:mi>
          <mml:msub>
            <mml:mi>&sigma;</mml:mi>
            <mml:mi>j</mml:mi>
          </mml:msub>
          </mml:msub>

  </mml:math>
 </inlineequation>

 is the background frequency of alphabet letter 
<inlineequation>  
<mml:math>
          <mml:msub>
            <mml:mi>&sigma;</mml:mi>
            <mml:mi>j</mml:mi>
          </mml:msub>

  </mml:math>
 </inlineequation>

. The normalizing
procedure has been left out in the above formula.
The purpose of using pre-determined background frequencies is to get a score which
depends on the column frequencies of the alphabet letters relative to a background distribution.
Note that the geometric mean interpretation does not hold for the symmetric version. The product
<inlineequation>  
<mml:math>

    <mml:mstyle displaystyle="true">
      <mml:munderover>
          <mml:mo>
&prod;    </mml:mo>
        <mml:mrow>
         <mml:mi>j</mml:mi>
          <mml:mo>=</mml:mo>
          <mml:mn>1</mml:mn>
        </mml:mrow>
        <mml:mi>Z</mml:mi>
      </mml:munderover>
      <mml:msup>

        <mml:mrow>
       <mml:mo     >(</mml:mo>
        <mml:mfrac>

        <mml:mrow>
     <mml:msubsup>
        <mml:mi>x</mml:mi>
        <mml:mi>k</mml:mi>
        <mml:mi>msa</mml:mi>
      </mml:msubsup>
      <mml:mo  stretchy="false">(</mml:mo>
      <mml:msub>
        <mml:mi>&sigma;</mml:mi>
        <mml:mi>j</mml:mi>
      </mml:msub>
      <mml:mo  stretchy="false">)</mml:mo>
        </mml:mrow>
         <mml:msub>
            <mml:mi>freq</mml:mi>
          <mml:msub>
            <mml:mi>&sigma;</mml:mi>
            <mml:mi>j</mml:mi>
          </mml:msub>
          </mml:msub>
        </mml:mfrac>
       <mml:mo  >)</mml:mo>

        </mml:mrow>
        <mml:mrow>
         <mml:msub>
            <mml:mi>e</mml:mi>
         <mml:msub>
            <mml:mi>&pi;</mml:mi>
            <mml:mi>i</mml:mi>

         </mml:msub>
          </mml:msub>
          <mml:mo stretchy="false">(</mml:mo>
          <mml:msub>
            <mml:mi>&sigma;</mml:mi>
            <mml:mi>j</mml:mi>
          </mml:msub>
          <mml:mo stretchy="false">)</mml:mo>
        </mml:mrow>
      </mml:msup>
    </mml:mstyle>
  </mml:math>
 </inlineequation>


cannot easily be given a probabilistic interpretation, because it contains an inverted view of the
sequences and the model, i.e. the emission probabilities are seen as a set of sequence residues, which are
emitted by the sequence frequency profile.

</para>
</sect3>
  <sect3 id="emission_scoring_method_dppi">
    <title>DPPI</title>

<para>

The fourth scoring method DPPI, is a combination of DP  ( <xref linkend="emission_scoring_method_dp"/> ) and PICASSO  ( <xref linkend="emission_scoring_method_picasso"/> ). Here, a sum of products between the two elements of
the two vectors
are calculated as in standard DP, but each term in the sum is also divided by the bakground frequency of that letter.





<equation>
  <mml:math>
    <mml:mstyle displaystyle="true">



    <mml:mrow>

      <mml:msubsup>
        <mml:mi>e</mml:mi>
        <mml:mi>q</mml:mi>
        <mml:mi>msa</mml:mi>
      </mml:msubsup>
      <mml:mo  stretchy="false">(</mml:mo>
      <mml:msubsup>
        <mml:mi>x</mml:mi>
        <mml:mi>k</mml:mi>
        <mml:mi>msa</mml:mi>
      </mml:msubsup>
      <mml:mo  stretchy="false">)</mml:mo>
      <mml:mo>=</mml:mo>


      <mml:munder>

          <mml:mo>
&sum;    </mml:mo>
        <mml:mi>j</mml:mi>

      </mml:munder>


<mml:mfrac><mml:mrow>
      <mml:msub>
        <mml:mi>e</mml:mi>
        <mml:mi>q</mml:mi>
      </mml:msub>
      <mml:mo  stretchy="false">(</mml:mo>
      <mml:msub>
        <mml:mi>&sigma;</mml:mi>
        <mml:mi>j</mml:mi>
      </mml:msub>
      <mml:mo  stretchy="false">)</mml:mo>
      <mml:mo>*</mml:mo>
      <mml:msubsup>
        <mml:mi>x</mml:mi>
        <mml:mi>k</mml:mi>
        <mml:mi>msa</mml:mi>
      </mml:msubsup>
      <mml:mo  stretchy="false">(</mml:mo>
      <mml:msub>
        <mml:mi>&sigma;</mml:mi>
        <mml:mi>j</mml:mi>
      </mml:msub>
      <mml:mo  stretchy="false">)</mml:mo>
</mml:mrow>

         <mml:msub>
            <mml:mi>freq</mml:mi>
          <mml:msub>
            <mml:mi>&sigma;</mml:mi>
            <mml:mi>j</mml:mi>
          </mml:msub>
          </mml:msub>

</mml:mfrac>



    </mml:mrow>
    </mml:mstyle>
  </mml:math>
</equation>





</para>
</sect3>
  <sect3 id="emission_scoring_method_smp">
    <title>SMP</title>


<para>
The fifth implementation of <inlineequation>  
<mml:math>
<mml:mrow>
 <mml:msup>
   <mml:mi>e</mml:mi>
   <mml:mi>msa</mml:mi>
 </mml:msup>
</mml:mrow>
  </mml:math>
 </inlineequation>
, SMP, is using a variant of the so called <emphasis>log-average score</emphasis>, which
originally is









<equation>
  <mml:math>
    <mml:mstyle displaystyle="true">



    <mml:mrow>

      <mml:msubsup>
        <mml:mi>e</mml:mi>
        <mml:mi>q</mml:mi>
        <mml:mi>msa</mml:mi>
      </mml:msubsup>
      <mml:mo  stretchy="false">(</mml:mo>
      <mml:msubsup>
        <mml:mi>x</mml:mi>
        <mml:mi>k</mml:mi>
        <mml:mi>msa</mml:mi>
      </mml:msubsup>
      <mml:mo  stretchy="false">)</mml:mo>
      <mml:mo>=</mml:mo>

      <mml:mo>log</mml:mo>


      <mml:munderover>
          <mml:mo>
&sum;    </mml:mo>
        <mml:mi>i=1</mml:mi>
<mml:mrow>
  <mml:mo stretchy="false">|</mml:mo>
          <mml:mo>
&sum;    </mml:mo>
  <mml:mo stretchy="false">|</mml:mo>
</mml:mrow>
      </mml:munderover>


      <mml:munderover>
          <mml:mo>
&sum;    </mml:mo>
        <mml:mi>j=1</mml:mi>
<mml:mrow>
  <mml:mo stretchy="false">|</mml:mo>
          <mml:mo>
&sum;    </mml:mo>
  <mml:mo stretchy="false">|</mml:mo>
</mml:mrow>
      </mml:munderover>

      <mml:msubsup>
        <mml:mi>x</mml:mi>
        <mml:mi>k</mml:mi>
        <mml:mi>msa</mml:mi>
      </mml:msubsup>
      <mml:mo  stretchy="false">(</mml:mo>
      <mml:msub>
        <mml:mi>&sigma;</mml:mi>
        <mml:mi>i</mml:mi>
      </mml:msub>
      <mml:mo  stretchy="false">)</mml:mo>

      <mml:mo>*</mml:mo>
      <mml:msub>
        <mml:mi>e</mml:mi>
        <mml:mi>q</mml:mi>
      </mml:msub>
      <mml:mo  stretchy="false">(</mml:mo>
      <mml:msub>
        <mml:mi>&sigma;</mml:mi>
        <mml:mi>j</mml:mi>
      </mml:msub>
      <mml:mo  stretchy="false">)</mml:mo>

<mml:mfrac><mml:mrow>

      <mml:msub>
        <mml:mi>p</mml:mi>
        <mml:mi>rel</mml:mi>
      </mml:msub>
      <mml:mo  stretchy="false">(</mml:mo>

        <mml:mi>i</mml:mi><mml:mo>,</mml:mo><mml:mi>j</mml:mi>

      <mml:mo  stretchy="false">)</mml:mo>





</mml:mrow>

<mml:mrow>
      <mml:msub>
        <mml:mi>p</mml:mi>
        <mml:mi>i</mml:mi>
      </mml:msub>


      <mml:msub>
        <mml:mi>p</mml:mi>
        <mml:mi>j</mml:mi>
      </mml:msub>
</mml:mrow>

</mml:mfrac>



    </mml:mrow>
    </mml:mstyle>
  </mml:math>
</equation>

where <inlineequation>  
<mml:math>

 <mml:msub>
   <mml:mi>p</mml:mi>
   <mml:mi>i</mml:mi>
 </mml:msub>

  </mml:math>
 </inlineequation>

 is the background distribution 
for letter  <inlineequation>  
<mml:math>
   <mml:mi>i</mml:mi>
 </mml:math>
 </inlineequation>

 and <inlineequation>  
<mml:math>

 <mml:msub>
   <mml:mi>p</mml:mi>
   <mml:mi>rel</mml:mi>
 </mml:msub>

  </mml:math>
 </inlineequation>

is a measure  of "relatedness" for letter pairs. The modhmm
implementation differs from this formula only in not taking the log of the double sum. A substitution matrix 
( <xref linkend="modhmm_subst_mtx_format"/> ) is supplied
by the user for the  <inlineequation>  
<mml:math>
<mml:mfrac><mml:mrow>

      <mml:msub>
        <mml:mi>p</mml:mi>
        <mml:mi>rel</mml:mi>
      </mml:msub>
      <mml:mo  stretchy="false">(</mml:mo>

        <mml:mi>i</mml:mi><mml:mo>,</mml:mo><mml:mi>j</mml:mi>

      <mml:mo  stretchy="false">)</mml:mo>





</mml:mrow>

<mml:mrow>
      <mml:msub>
        <mml:mi>p</mml:mi>
        <mml:mi>i</mml:mi>
      </mml:msub>


      <mml:msub>
        <mml:mi>p</mml:mi>
        <mml:mi>j</mml:mi>
      </mml:msub>
</mml:mrow>

</mml:mfrac>


  </mml:math>
 </inlineequation>
 part.
This scoring method increases the time complexity of the forward and backward algorithms with a factor the size of the
alphabet to the power of 2.

</para>
</sect3>
  <sect3 id="emission_scoring_method_smdp">
    <title>SMDP</title>


<para>
The sixth implementation of <inlineequation>  
<mml:math>
<mml:mrow>
 <mml:msup>
   <mml:mi>e</mml:mi>
   <mml:mi>msa</mml:mi>
 </mml:msup>
</mml:mrow>
  </mml:math>
 </inlineequation>, SMDP, is a faster version of SMP  ( <xref linkend="emission_scoring_method_smp"/> ), which puts more focus on the query sequence.






<equation>
  <mml:math>
    <mml:mstyle displaystyle="true">



    <mml:mrow>

      <mml:msubsup>
        <mml:mi>e</mml:mi>
        <mml:mi>q</mml:mi>
        <mml:mi>msa</mml:mi>
      </mml:msubsup>
      <mml:mo  stretchy="false">(</mml:mo>
      <mml:msubsup>
        <mml:mi>x</mml:mi>
        <mml:mi>k</mml:mi>
        <mml:mi>msa</mml:mi>
      </mml:msubsup>
      <mml:mo  stretchy="false">)</mml:mo>
      <mml:mo>=</mml:mo>


      <mml:munderover>
          <mml:mo>
&sum;    </mml:mo>
        <mml:mi>j=1</mml:mi>
<mml:mrow>
  <mml:mo stretchy="false">|</mml:mo>
          <mml:mo>
&sum;    </mml:mo>
  <mml:mo stretchy="false">|</mml:mo>
</mml:mrow>
      </mml:munderover>

      <mml:msubsup>
        <mml:mi>x</mml:mi>
        <mml:mi>k</mml:mi>
        <mml:mi>msa</mml:mi>
      </mml:msubsup>
      <mml:mo  stretchy="false">(</mml:mo>
      <mml:msub>
        <mml:mi>&sigma;</mml:mi>
        <mml:mi>query</mml:mi>
      </mml:msub>
      <mml:mo  stretchy="false">)</mml:mo>

      <mml:mo>*</mml:mo>
      <mml:msub>
        <mml:mi>e</mml:mi>
        <mml:mi>q</mml:mi>
      </mml:msub>
      <mml:mo  stretchy="false">(</mml:mo>
      <mml:msub>
        <mml:mi>&sigma;</mml:mi>
        <mml:mi>j</mml:mi>
      </mml:msub>
      <mml:mo  stretchy="false">)</mml:mo>

<mml:mfrac><mml:mrow>

      <mml:msub>
        <mml:mi>p</mml:mi>
        <mml:mi>rel</mml:mi>
      </mml:msub>
      <mml:mo  stretchy="false">(</mml:mo>        <mml:mi>query</mml:mi><mml:mo>,</mml:mo><mml:mi>j</mml:mi>

      <mml:mo  stretchy="false">)</mml:mo>





</mml:mrow>

<mml:mrow>
      <mml:msub>
        <mml:mi>p</mml:mi>
        <mml:mi>query</mml:mi>
      </mml:msub>


      <mml:msub>
        <mml:mi>p</mml:mi>
        <mml:mi>j</mml:mi>
      </mml:msub>
</mml:mrow>

</mml:mfrac>



    </mml:mrow>
    </mml:mstyle>
  </mml:math>
</equation>


where 


<inlineequation>  
<mml:math>


      <mml:msubsup>
        <mml:mi>x</mml:mi>
        <mml:mi>k</mml:mi>
        <mml:mi>msa</mml:mi>
      </mml:msubsup>
      <mml:mo  stretchy="false">(</mml:mo>
      <mml:msub>
        <mml:mi>&sigma;</mml:mi>
        <mml:mi>query</mml:mi>
      </mml:msub>
      <mml:mo  stretchy="false">)</mml:mo>


  </mml:math>
 </inlineequation>
 represents the alignment share of the query sequence's letter at
position/column 

<inlineequation>  
<mml:math>
        <mml:mi>k</mml:mi>
  </mml:math>
 </inlineequation>


. This scoring method is an approximation of the SMP method which has the advantage of
being a factor the size of the alphabet faster.
</para>
</sect3>
  <sect3 id="emission_scoring_method_smdpp">
    <title>SMDPP</title>
<para>
The seventh implementation, SMDPP, is essentially the same as the SMDP method ( <xref linkend="emission_scoring_method_smdp"/> ), but for all letters in a
certain position for which the substitution matrix value is below a certain threshold, the sequence column
shares are evened out, i.e. redistributed uniformly over these letters. The objective is to disregard
the significance of a random, "unlikely" mutation.




    </para>
  </sect3>
  </sect2>
  </sect1>


<bibliography>
    <title>References</title>

<!--
    <biblioentry id="bib.xml">
    <abbrev id="bib.xml.abbrev">XML98</abbrev>
   <authorgroup>
 <author><firstname>A.</firstname><surname>Heger</surname></author>
 <author><firstname>L.</firstname><surname>Holm</surname></author>
    </authorgroup>
    <title>Exhaustive enumeration of protein domain families</title>
    <publisher>
 <publishername>J. Mol. Biol.</publishername>
    </publisher>
   <year>2003</year>
   <volumenum>328</volumenum>
 <pagenums>749-767</pagenums>
    </biblioentry>
-->
<bibliomixed id="heger:2003">


  <bibliomset relation="article">
 

  

    <author><firstname>A.</firstname> <surname>Heger</surname></author> and <author><firstname>L.</firstname> <surname>Holm</surname></author>
    <title role="article">Exhaustive enumeration of protein domain families</title>.
  </bibliomset>
  <bibliomset relation="journal">
    <title>J. Mol. Biol.</title> 
    <volumenum>328</volumenum>:<pagenums>749-767</pagenums>, <pubdate>2003</pubdate></bibliomset>.
</bibliomixed>


<bibliomixed  id="hughey:1996">
  <bibliomset relation="article">
 
 <abbrev id="hughey:1996.abbrev">Hughey:1994</abbrev>


    <author><firstname>R.</firstname> <surname>Hughey</surname></author> and <author><firstname>A.</firstname> <surname>Krogh</surname></author>
    <title role="article">Hidden Markov models for sequence analysis: extension and analysis of the basic method</title>.
  </bibliomset>
  <bibliomset relation="journal">
    <title>Comput. Appl. Biosci.</title>.,
    <volumenum>12</volumenum>(<issuenum>2</issuenum>):<pagenums>95-107</pagenums>, <pubdate>1996</pubdate></bibliomset>.
</bibliomixed>





<bibliomixed  id="krogh:1994">
  <bibliomset relation="article">
 


<author><firstname>A.</firstname> <surname>Krogh</surname></author>
    <title role="article">Hidden Markov models for labeled sequences</title>.
  </bibliomset>
  <bibliomset relation="journal">
    <title> Proceedings of the 12th IAPR International Conference on
  Pattern Recognition</title>., Los Alamitos, California,
   pages <pagenums>140-144</pagenums>, <pubdate>1994</pubdate></bibliomset> <publisher><publishername>IEEE  Computer Society Press</publishername></publisher>.
</bibliomixed>









<bibliomixed id="mittelman:2003">
  <bibliomset relation="article" >
 


    <author><firstname>D.</firstname> <surname>Mittelman</surname></author>, <author><firstname>R.</firstname> <surname>Sadreyev</surname></author> and  <author><firstname>N.</firstname> <surname>Grishin</surname></author>
    <title role="article">Probabilistic scoring measures for profile-profile comparison yield
  more accurate short seed alignments</title>.
  </bibliomset>
  <bibliomset relation="journal">
    <title>Bioinformatics</title>.,
    <volumenum>19</volumenum>:<pagenums>1531-1539</pagenums>, <pubdate>2003</pubdate></bibliomset>.
</bibliomixed>




<bibliomixed id="sjolander:1996">
  <bibliomset relation="article" >
 



    <author><firstname>K.</firstname> <surname>Sjölander</surname></author>, <author><firstname>K.</firstname> <surname>Karplus</surname></author>, <author><firstname>R.</firstname> <surname>Hughey</surname></author>,     <author><firstname>A.</firstname> <surname>Krogh</surname></author>,  <author><firstname>I.S.</firstname> <surname>Mian</surname></author> and  <author><firstname>D.</firstname> <surname>Haussler</surname></author>
    <title role="article">Dirichlet mixtures: a method for improved detection of weak but
  significant protein sequence homology</title>.
  </bibliomset>
  <bibliomset relation="journal">
    <title>Comput. Appl. Biosci.</title>.,




    <volumenum>12</volumenum>(<issuenum>4</issuenum>):<pagenums>327-345</pagenums>, <pubdate>1996</pubdate></bibliomset>.


</bibliomixed>

    </bibliography>


</article>




