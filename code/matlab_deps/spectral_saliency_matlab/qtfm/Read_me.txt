Quaternion toolbox for Matlab - installation instructions
---------------------------------------------------------

To install this toolbox:

1. Unzip the distribution file and move the directory/folder
   to a convenient location. This does not need to be in the
   same location as toolboxes supplied with Matlab, although
   it could be if your Matlab installation is to be used by
   multiple users on the same machine.

2. Set the Matlab path to include the directory/folder qtfm.
   This is done from the Matlab File -> Set Path menu.
   [The QTFM folder must be near the top of the path, i.e.
   higher than the standard Matlab folders, otherwise the
   overloading of Matlab functions will not work.]

3. Help information is available from the Matlab Start menu
   (at bottom left of the Matlab window) under Toolboxes
   -> Quaternion; or in the Matlab Help browser (Quaternion
   Toolbox should appear as one of the help items in the left
   side pane). You can also access help information in the
   Matlab command window: try help qtfm.

4. help quaternion/xxx shows help text for each function xxx
   and more detailed help is available in the help documentation
   (see 3 above).

5. To run a test of the toolbox, set the current directory to
   .../qftm/test and type the command 'test' in the Matlab
   command window. This runs a test of many parts of the
   toolbox and will allow you to confirm that installation is
   correct. It is also possible to invoke the test code from
   the Matlab start menu (Toolboxes -> Quaternion -> Run test
   code).

6. If you find the toolbox useful and would like to be kept
   informed of updates and future releases, subscribe to the
   mailing list qtfm-announce@lists.sourceforge.net. You can
   subscribe to this list at:
   https://lists.sourceforge.net/lists/listinfo/qtfm-announce

7. LaTeX/BiBTeX users please note the file qtfm.bib provided
   in the directory tex/bibtex if you wish to cite QTFM itself
   or one of the published papers cited within the source code.

8. We would welcome contributed code, ideas for new functions,
   or expanded examples to include in the help documentation.
   (Contributed code must be made available under the same
   license terms as QTFM itself). Please contact us by email
   as in below.

9. To send feedback, code, documentation corrections or bug
   reports, please use the reporting mechanisms at the project
   webpage at Sourceforge

   https://sourceforge.net/projects/qtfm/

   or send email to:

   sangwine@users.sourceforge.net
   n-le_bihan@users.sourceforge.net

Steve Sangwine
Nicolas Le Bihan

27 March 2006
Updated 5 June 2006 and 23 May 2008.
$Id: Read_me.txt,v 1.6 2009/12/24 10:34:25 sangwine Exp $
