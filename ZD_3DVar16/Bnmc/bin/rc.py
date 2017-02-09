#
# Python module 'rc' .
#

# ------------------------------------------------
# help
# ------------------------------------------------

"""
NAME
    rc - python module to read settings from a file with X-Resources syntax

RCFILES

    A rcfile is a text file with key/value pairs seperated by a ':', e.g.

        my.flag    :  T
        my.answer  :  42

    The following functionality is supported:

    * Line continuation
     
        Long values could be continued at the next line after a '\' as last character.
     
    * Annotation
     
        o Empty lines are ignored.

        o Comment lines are introduced by a '!' as first character.

        o Comment starting with '!' is stripped from the values.
          To have a value including exclamation marks, use '\!' but do
          not expect that the rest of the value is scanned for comment too:
       
            my.value      :   -999    ! just an integer value
            my.message    :   This value has 64 characters \! Count if you don't believe it ...

    * Variable subsitituion
     
        Variables are enclosed by '${' and '}', e.g. '${HOME}' .
        The python functions also support other patterns,
        which should then be specified with the optional 'marks' tupple.

        In case some substitutions could not be revold an error message will be raised.             

        o Substitute of environment variables:

            tm.dir : ${HOME}/TM5/cy3

        o Substitute values of other keys in a line:

            build.dir            :  ${tm.dir}/build
            grid                 :  glb300x200
            input.${grid}.path   :  /data/input/${grid}

          Substitions are allowed in both key names as well as values.
          The substitutions are performed in a loop until nothing
          has to be substituted anymore, or some substitutions could
          not be applied at al; for the later an error is raised.
          Values to be substituted could therefore be set before and
          after they are used.

          Note that if a key has the same name as an environment variable,
          the new value will be assigned to the key instead of the value
          retrieved from the environment:

            HOME      :  /some/other/dir/

        o Substitude some specials:

            ${pid}     # evaluates to the current process id; 
                       # useful for names of log files etc
            ${script}  # evaluates to the base name of the calling script, 
                       # thus without .py etc

    * Evaluation of expressions is applied to all values enclosed by '$((..))' .
      The enclosed value should be a valid Python expression:
      
        ntask       :  4
        nthread     :  2
        ncore       :  $(( ${ntask} * ${nthread} ))

    * Special processing
     
        o Include the key/value pairs from another file:

            #include an/other.rc

        o Conditional expressions:

            #if ${my.number} == 1
            message    : Welcome
            #else
            message    : Whatever ...
            #endif

          The conditions should be valid python expressions that evaluate to a boolean;
          value substitutions are performed before evaluation. Examples:

            ${my.runmode} == 4
            "${my.tracer}" == "CH4"

          Keep it simple! Very complicated and nested if-statements might not be
          resolved correctly, and are in any case not easy to understand for other users!

          In the example above, an exception could be raised by the special error expression;
          everything behind the '#error' mark is displayed as an error message:

            #error No settings provided for value : ${my.value}
              
        o A for-loop could be used to quickly set a number of similar settings:
        
            #for XX in AA BB CC :
            setting.XX   :  This is the value for XX.
            #endfor
        
          will expand to:
          
            setting.AA   :  This is the value for AA.
            setting.BB   :  This is the value for BB.
            setting.CC   :  This is the value for CC.
       

USAGE AS PYTHON MODULE

    Import the module with:
    
        import rc

    Initialiase by reading all settings in a rcfile,
    supporting the functionality described in the 'RCFILES' section.

        rcf = RcFile( 'settings.rc' )

    The initialisation accepts some optional arguments.
    Set the silent flag to True to ignore warnings.

        rcf = RcFile( 'settings.rc', silent=False )

    Use the optional 'marks' tupple to define that variables to be expanded
    are marked other than '${..}' but rather '<mark1>..<mark2>' :

        rcf = RcFile( 'settings.rc', marks=('${','}') )

    Test to see if a key is defined:

        if rcf.has_key('my.flag') :
            print 'value of my flag is : ', rcf['my.flag']
    
    Extract a list with all keys:
    
        rcf.keys()

    A 'get' function is provided to extract values:

     * by default, the 'get' function returns the value as a str type:

         s = rcf.get('my.value')
    
     * a second argument is the name of the python type to which
       the value is converted to:
    
         i = rcf.get('my.flag','int')

     * if the return value should be a 'bool', the result is
         True  for values     : 'True' , 'T', 'yes', or '1' ,
         and False for value  : 'False', 'F', 'no' , or '0' ;
       for other values an error is raised;
     
     * return a default value if the key is not found:

            rcf.get( 'my.flag', default=False )
    
     * print a debug message to the logging system for each extracted key:
     
            rcf.get( 'my.flag', verbose=True ) 

    Add a new value, comment is optional:

        rcf.add( 'my.iter', 2, comment='iteration number for restart' )

    Assign a new value to an existing key:

        rcf.replace( 'my.flag', True )

    Scan a character line for all occurances of ${<key>} and subsitute for
    the rc value assigned to <key> :

        line = rcf.substitute( line )

    Write the dictionary (with all variables expanded and included files included)
    to new file:

         rcf.WriteFile('newfile.rc')

         
USAGE AS PYTHON MODULE - BACKWARDS COMPATIBILITY

    For backwards compatibility with older implementations of the rc.py module,
    two extra routines are available.
    
    To read rc-file by making an instance of the RcFile class, 
    and to returns a dictionary of values only, use:
            
        rcdict = read( 'test.rc' [,silent=False] )
        
    Create a new rcfile and fill with key/values from a dictionary:

        write( 'test.rc', rcdict )


USAGE VIA SCRIPT

    Use the 'rcget' script to extract values from a rcfile in
    non-python application.


HISTORY

    2008? Andy Jacobson, NOAA
      Translation to python of original shell script 'go_readrc' .
    2009-06 Wouter Peters, WUR
      Support substitution of previously defined variables.
    2009-06 Arjo Segers, TNO
      Support include files.
    2009-09 Arjo Segers, TNO
      Re-coded into class.
      Implemented substitution loop.
    2009-11 Arjo Segers, JRC
      Added main program to run this file as a shell script.
      Added replace and substitute routines.
    2010-03 Arjo Segers, JRC
      Support simple if-statements.
      Support comment in values.
    2010-07 Wouter Peters, WUR
      Downgraded to work for python 2.4.3 too.
      Added read/write routines for backwards compatibility.
    2010-07-27 Arjo Segers, JRC
      Maintain list with rcfile names and line numbers to be displayed
      with error messages to identify where problematic lines are found.
    2010-07-28 Andy Jacobson, NOAA
      Add second dictionary of key,linetrace values to help track the 
      provenance of #included keys (to debug multiple key instances).
      Identify duplicate keys by checking on different source lines
      instead of checking if the values are different.
    2010-10 Arjo Segers
      Restructured processing using classes for source lines
      and rcfile values, and resolve using recursive calls.
      Added evaluation of expression enclosed by $((.)) .
      Added for-loop.
      Removed main program and stored this in the auxilary script 'rcget' .
"""


# ------------------------------------------------
# logging
# ------------------------------------------------

# modules:
import sys
import logging

# setup logging if not done yet:
logging.basicConfig( stream=sys.stdout, level=logging.INFO,
                       format='[%(levelname)-8s] %(message)s' )

# get root logger instance:
logger = logging.getLogger()


# ------------------------------------------------
# classes
# ------------------------------------------------


class RcFile( object ) :

    """
    Class to store settings read from a rcfile.
    """

    def __init__( self, filename, raw=False, silent=False, marks=('${','}'), debug=False ) :

        """ 
        
        Usage:
        
          rcf = RcFile( 'settings.rc' [,raw=False], [,silent=False] [marks=('${','}')] )

        Read rc-file and store content.
        Variable substitutions are applied and special lines are evaluated,
        unless 'raw' is set to True.
        Key/value pairs are stored and can be retrieved using the 'get' function.
        Do not shout messages if silent is set to True. 
        The 2-item tupple (mark1,mark2) could be used to re-define the default
        substitution pattern '${..}' into something else:
          <mark1>...<mark2>

        """

        # external:
        import re
        import os
        import sys
        import logging
        
        # info ...
        if debug :
            logger.debug( 'reading rcfile %s ...' % filename )
        #endif

        # check ...
        if not os.path.exists(filename) :
            msg = 'rcfile not found : %s' % filename ; logger.error(msg)
            raise IOError, msg
        #endif
        
        # store file name:
        self.filename = filename
        # store rc-file root directory:
        self.rootdir = os.path.split(filename)[0]
        
        # read lines from file:
        SL = RcLines()
        SL.ReadFile( filename, raw=raw )

        # collect into tree:
        self.tree = RcLines()
        self.tree.ReadTree( SL )

        # show initial listing:
        if debug :
            logger.info( '---[rcfile]------------------------------------------' )
            self.tree.Show( level=logging.INFO, indent='', autoindent=False )
            logger.info( '-----------------------------------------------------' )
            logger.info( '' )
        #endif

        # raw mode ?
        if raw :
        
            # extract current dictionary:
            self.rcvalues = self.tree.Extract()
            # nothing undefined:
            undefined_keys = []
            
            # show if needed:
            if debug :
                logger.info( '---[rcfile after raw]--------------------------------' )
                self.tree.Show( level=logging.INFO, indent='', autoindent=False )
                logger.info( '-----------------------------------------------------' )
                logger.info( '' )
            #endif
            
        else :

            # start with empty rc dictionary:
            self.rcvalues = {}
            # initial pass conter:
            ipass = 0
            # until all resolved:
            while True :

                # increase counter:
                ipass = ipass + 1

                # reset changed flag:
                self.tree.setflag( changed=False )

                # try to resolve:
                undefined_keys = self.tree.resolve( self.rcvalues )

                # try to evaluate:
                self.tree.evaluate()

                # show current listing:
                if debug :
                    logger.info( 'pass %i ...' % ipass )
                    logger.info( '---[rcfile]------------------------------------------' )
                    self.tree.Show( level=logging.INFO, indent='', autoindent=False )
                    logger.info( '-----------------------------------------------------' )
                    logger.info( '' )
                #endif

                # extract current dictionary:
                self.rcvalues = self.tree.Extract()

                # no changes during latest 'resolve' ? then leave:
                if (ipass > 1) and (not self.tree.changed()) :
                    if debug : logger.info( 'nothing changed during last resolve; break ...' )
                    break
                #endif

                # for safety ...
                if ipass == 100 :
                    logger.error( 'resolving rc file has reached pass %i ; something wrong ?' % ipass )
                    break
                #endif

                # testing ...
                #break
            #endwhile
            
        #endif  # raw

        # not completely resolved yet ?
        if not self.tree.resolved() :
            logger.error( 'could not resolve rcfile ...' )
            # show unresolved special lines:
            logger.error( '  unresolved lines:' )
            self.tree.Show( level=logging.ERROR, indent='    ', unresolved=True )
            # any undefined keys ?
            if len(undefined_keys) > 0 :
                # list with rcvalues that could not be resolved due to
                # unresolved substitutions:
                kwavs = self.tree.KeysWithUnresolvedValue()
                # info :
                logger.error( '  undefined variables : ' )
                # loop over undefined keys:
                for key in undefined_keys :
                    # in list with keys that could not be resolved ? then skip:
                    if key in kwavs : continue
                    # display:
                    logger.error( '    %s' % key )
                    self.tree.Show( indent='      ', undefined_key=key )
                #endfor
            #endif
            # quit:
            raise Exception
        #endif

        # show dictionary:
        if debug :
            logger.info( '' )
            logger.info( '---[rcvalues]----------------------------------------' )
            for key in self.rcvalues.keys() :
              logger.info( '%s = %s' % (key,self.rcvalues[key].getvalue()) )
            #endfor
            logger.info( '-----------------------------------------------------' )
            logger.info( '' )
        #endif

    #enddef  # __init__
    
    
    # ***
    
    
    def has_key( self, key ) :
    
        # from dictionary:
        return self.rcvalues.has_key(key)
        
    #enddef
    
    
    # ***
    
    
    def keys( self ) :
    
        # from rcvalues dictionary:
        return self.rcvalues.keys()
        
    #enddef
    
    
    # ***


    def get( self, key, totype='', default=None, verbose=False ) :
    
        """
        rcf.get( 'my.value' [,default=None] )
        Return element 'key' from the dictionary.
        If the element is not present but a default is specified, than return
        the default value.
        If 'verbose' is set to True, then print debug messages to the logging
        about which values is returned for the given key.
        The option argument 'totype' defines the conversion to a Python type.
        If 'totype' is set to 'bool', the return value is the
        boolean True for values 'T', 'True', 'yes', and '1',
        and False for 'F', 'False', 'no', or '0' ;
        for other values, an exception will be raised.
        """
        
        # element found ?
        if self.rcvalues.has_key(key) :
            # copy value:
            value = self.rcvalues[key].getvalue()
            # convert ?
            if totype == 'str' :
                # leave as it is ...
                pass
            elif totype == 'bool' :
                # convert to boolean:
                if value in ['T','True','yes','1'] :
                    value = True
                elif value in ['F','False','no','0'] :
                    value = False
                else :
                    logger.error( "value of key '%s' is not a boolean : %s" % (key,str(value)) )
                    raise Exception
                #endif
            elif len(totype) > 0 :
                # convert to other type ...
                value = eval( '%s(%s)' % (totype,value) )
            #endif
            # common mistake ...
            if (value in ['True','False','T','F']) and (len(totype)==0) :
                logger.error( 'Value for "%s" is "%s", but not explicitly asked for a bool ...' % (key,value) )
                raise Exception
            #endif
            # for debugging ...
            if verbose : logger.debug( 'rc setting "%s" : "%s"' % (key,str(value)) )
        else :
            # default value specified ?
            if default != None :
                # copy default:
                value = default
                # for debugging ...
                if verbose : logger.debug( 'rc setting "%s" : "%s" (deault)' % (key,str(value)) )
            else :
                # something wrong ...
                logger.error( "key '%s' not found in '%s' and no default specified" % (key,self.filename) )
                raise Exception
            #endif
        #endif
        
        # ok
        return value
        
    #enddef
    
    
    # ***
    
    
    def replace( self, key, val ) :
    
        """
        Replace the value of a key by a new value.
        """
        
        # check ...
        if not self.rcvalues.has_key(key) :
            logger.error( 'could not replace value, key not found : %s' % key )
            raise ValueError
        #endif
        # apply replace to dictionairy with values:
        self.rcvalues[key].replacevalue( val )
        # also the key:value line the tree with source lines:
        self.tree.replacevalue( key, val )
    
    #enddef
    
    
    # ***
    
    
    def add( self, key, val, comment='' ) :
    
        """Add a new key/value pair."""
        
        # add lines; use empty filename and zero line numbers:
        # ~ empty line:
        rcline = RcLine( '', '', 0 )
        self.tree.append( rcline )
        # ~ comment ?
        if len(comment) > 0 :
            rcline = RcLine( '! %s' % comment, '', 0 )
            self.tree.append( rcline )
        #endif
        # ~ rcline with key:value pair:
        rcline = RcLine( '%s : %s\n' % (key,str(val)), '', 0 )
        self.tree.append( rcline )

        # add to dictionary:
        self.rcvalues.append( RcValue( key, val, rcline ) )
        
        # ok
        return
    
    #enddef
    
    
    # ***
    
    
    def substitute( self, line, marks=('${','}') ) :
    
        """
        Return a line with all '${..}' parts replaced by the corresponding rcfile values.
        The 2-item tupple (mark1,mark2) could be used to re-define the default
        key pattern '${..}' into something else:
          <mark1>...<mark2>
        """
        
        # external:
        import re
        
        # ensure that common marks are evaluated correctly:
        start_mark = marks[0].replace('{','\{').replace('<','\<').replace('$','\$')
        close_mark = marks[1].replace('}','\}').replace('>','\>')

        # set syntax of keywords to be matched, e.g. '${...}' :
        pattern = start_mark+'[A-Za-z0-9_.-]+'+close_mark

        # make a regular expression that matches all variables:
        rc_varpat = re.compile( pattern )

        # search all matching paterns:
        pats = re.findall(rc_varpat,line)
        # loop over matches:
        for pat in pats :
            # remove enclosing characters:
            key = pat.lstrip(start_mark).rstrip(close_mark)
            # test dictionary for matching key:
            if self.rcvalues.has_key(key) :
                # get previously defined value:
                val = self.rcvalues[key].getvalue()
                # substitute value:
                line = line.replace(pat,val)
            #endif
        #endfor  # matched patterns

        # ok
        return line
        
    #enddef
    

    # ***


    def getlines( self ) :

        """
        Return list with processed rcfile lines.
        """
        # ok
        return self.tree.getlines()
    #endif
    

    # ***


    def WriteFile( self, filename ) :

        """ write the dictionary to file"""

        # get lines:
        lines = self.tree.getlines()
        # open file for writing:
        f = open(filename,'w')
        # write line plus newline:
        for line in lines : f.write( line+'\n' )
        # close file:
        f.close()
        
    #endif
    

#endclass    # RcFile


#-------------------------------------------------
# rcline(s)
#-------------------------------------------------


class RcLine( object ) :

    """
    Content and traceback info on a source line.
    Data fields:
        line        : content
        filename    : name of the file with the line
        lineno      : line number
        lineno2     : extra line number in case of continued lines
    """
    
    # NOTE: class(str) does not work, therefore .line attribute
    
    def __init__( self, line, filename, lineno, lineno2=None ) :
        """
        Define a source line and traceback info.
        """
        # new string:
        self.line = line
        # store traceback info:
        self.filename = filename
        self.lineno = lineno
        self.lineno2 = lineno
        if lineno2 != None : self.lineno2 = lineno2
        # comment or white line ?
        self._comment = (len(self.line.strip()) == 0) or \
                        self.line.strip().startswith('!')
        # no need to resolve comment lines:
        self._resolved = self._comment
        # init changed flag:
        self._changed  = False
        # special line?
        self._special = self.line.startswith('#')
        self._hold    = False   # wait with resolving until if-condition is evaluated etc
    #enddef   # __init__
    
    # *
    
    def extend( self, line, lineno ) :
        """
        Add two source lines together.
        Lines are connected with a single white space.
        """
        # add lines:
        self.line = self.line.strip()+' '+line.strip()
        # store second line number:
        self.lineno2 = lineno
    #enddef
    
    # *
    
    def startswith(self,*args) :
        return self.line.startswith(*args)
    #endif
    
    # *
    
    def strip(self,*args) :
        return self.line.strip(*args)
    #endif
    
    # *
    
    def replace(self,*args) :
        """
        Replace pattern in source line.
        """
        self.line = self.line.replace(*args)
    #endif
    
    # *
    
    def replacevalue(self,key,value) :
        """
        Replace original value assigned to key by new value.
        """
        # is this actually a 'key:value' line ?
        if ':' in self.line :
            # split in key and value
            okey,ovalue = self.line.split(':',1)
            # match ?
            if okey.strip() == key :
                # leave original key, replace value:
                self.line = '%s: %s' % (okey,value)
            #endif
        #endif
    #endif
    
    # *
    
    def copy( self ) :
        """
        Return a copy.
        """
        # create new instance:
        res = RcLine( self.line, self.filename, 
                            self.lineno, lineno2=self.lineno2 )
        # copy flags:
        res._comment  = self._comment
        res._resolved = self._resolved
        res._changed  = self._changed
        res._special  = self._special
        res._hold     = self._hold
        # ok
        return res
    #enddef
    
    # *
    
    def getline( self ) :
        """
        Return source line.
        """
        return self.line
    #enddef

    # *
    
    def setline( self, line ) :
        """
        Fill source line.
        """
        self.line = line
    #enddef

    # *
    
    def __str__( self ) :
        """
        Format sourceline and traceback info.
        """
        # single or multiple lines ?
        if self.lineno2 == self.lineno :
            # display as :  "test.src"     123 | This is a line
            msg = '"%s" %7i %s | %s' % (self.filename,self.lineno,self.flags(),str.__str__(self.line))
        else :
            # display as :  "test.src" 123-125 | This is a line
            msg = '"%s" %3i-%3i %s | %s' % (self.filename,self.lineno,self.lineno2,self.flags(),str.__str__(self.line))
        #endif
        # ok
        return msg
    #enddef  # __str__
        
    # *
    
    def traceback( self ) :
        """
        Format traceback info.
        """
        # single or multiple lines ?
        if self.lineno2 == self.lineno :
            # display as :  line 1234 in "test.src"
            msg = 'line %i in "%s"' % (self.lineno,self.filename)
        else :
            # display as :  lines 1234-1235 in "test.src"
            msg = 'lines %i-%i in "%s"' % (self.lineno,self.lineno2,self.filename)
        #endif
        # ok
        return msg
    #enddef  # traceback
    
    # *
    
    def flags( self ) :
        """
        Return character list to denote if some internal flag is set or not;
        a character is either '-' or (in this order):
         'c' comment, 's' special, 'h' holded, 'r' resolved, '*' changed .
        """
        res = ''
        if self._comment : res = res+'c'
        else : res = res+'-'
        if self._special : res = res+'s'
        else : res = res+'-'
        if self._hold : res = res+'h'
        else : res = res+'-'
        if self._resolved : res = res+'r'
        else : res = res+'-'
        if self._changed : res = res+'*'
        else : res = res+'-'
        return res
    #endif
        
    # *
    
    def setflag( self, empty=None, hold=None, changed=None, resolved=None ) :
        """
        Set internal flag.
        """
        if empty != None :
            self._comment = empty
            if self._comment :
                self._special  = False
                self._hold     = False
                self._resolved = True
            #endif
        #endif
        if hold     != None : self._hold     = hold
        if changed  != None : self._changed  = changed
        if resolved != None : self._resolved = resolved
    #enddef
    
    # *
    
    def commented( self ) :
        """
        Return True if this is an empty line (white line or comment).
        """
        # copy value:
        return self._comment
    #enddef

    # *
    
    def resolved( self ) :
        """
        Return True if line is resolved.
        """
        # copy value:
        return self._resolved
    #enddef
    
    # *
    
    def special( self ) :
        """
        Return True if this is a special line (starts with #).
        """
        # copy value:
        return self._special
    #enddef

    # *
    
    def changed( self ) :
        """
        Return True if some changes were applied 
        during the latest call to 'resolve()'.
        """
        # copy value:
        return self._changed
    #enddef
    
    # *
    
    def resolve( self, rcd, marks=('${','}') ) :

        """
        Resolve line using key/values from dictionary 'rcd'.
        Return list with undefined keys.
        """
        
        # modules:
        import sys
        import os
        import re
        
        # init result:
        undefined_keys = []
        
        # not hold and not resolved yet ?
        if (not self._comment) and (not self._hold) and (not self._resolved) :
        
            # ensure that common marks are evaluated correctly:
            start_mark = marks[0].replace('{','\{').replace('<','\<').replace('$','\$')
            close_mark = marks[1].replace('}','\}').replace('>','\>')

            # set syntax of keywords to be matched, e.g. '${...}' :
            pattern = start_mark+'[A-Za-z0-9_.-]+'+close_mark

            # make a regular expression that matches all variables:
            rc_varpat = re.compile( pattern )

            # search all matching paterns:
            pats = re.findall(rc_varpat,self.line)
            # counter for unexpanded substitutions:
            ntobedone = 0
            # loop over matches:
            for pat in pats :
                # remove enclosing characters:
                key = pat.lstrip(start_mark).rstrip(close_mark)
                # test some dictionaries for matching key:
                if rcd.has_key(key) :
                    # get previously defined value:
                    val = rcd[key].value
                    # substitute value:
                    self.line = self.line.replace(pat,val)
                    # something was done:
                    self.setflag(changed=True)
                elif os.environ.has_key(key) :
                    # get value from environment:
                    val = os.environ[key]
                    # substitute value:
                    self.line = self.line.replace(pat,val)
                    # something was done:
                    self.setflag(changed=True)
                elif key == 'pid' :
                    # special value: process id; convert to character:
                    val = '%i' % os.getpid()
                    # substitute value:
                    self.line = self.line.replace(pat,val)
                    # something was done:
                    self.setflag(changed=True)
                elif key == 'script' :
                    # special value: base name of the calling script, without extension:
                    script,ext = os.path.splitext(os.path.basename(sys.argv[0]))
                    # substitute value:
                    self.line = self.line.replace(pat,script)
                    # something was done:
                    self.setflag(changed=True)
                else :
                    # could not substitute yet; set flag:
                    ntobedone = ntobedone + 1
                    # add to list with unresolved keys:
                    if key not in undefined_keys : undefined_keys.append(key)
                #endif
            #endfor  # matched patterns

            # resolved if all substituted, unless special line:
            self._resolved = ntobedone == 0
            
            # if status changed to resolved, than that is a change:
            if self._resolved : self.setflag(changed=True)
                
        #endif  # not resolved yet
        
        # ok:
        return undefined_keys
        
    #enddef    # resolve

    # *
    
    def evaluate( self ) :
        """
        Evaluate $((..)) expresions.
        """

        # modules:
        import os
        import re
        
        # not hold, but all substitutions resolved ?
        if (not self._comment) and (not self._hold) and self._resolved :

            # set syntax of keywords to be matched, e.g. '${...}' :
            pattern = '\$\(\(.*\)\)'

            # make a regular expression that matches all variables:
            rc_varpat = re.compile( pattern )

            # search all matching paterns:
            pats = re.findall(rc_varpat,self.line)
            # loop over matches:
            for pat in pats :
                # remove enclosing characters:
                expression = pat.lstrip('$((').rstrip('))')
                # evaluate:
                value = eval( expression )
                # substitute value:
                self.line = self.line.replace(pat,str(value))
                # this is a change ...
                self.setflag( changed=True )
            #endfor  # matched patterns
                
        #endif  # expressions could be evaluated
        
    #enddef
    
    # *
    
    def comment( self ) :
        """
        Comment the line.
        """
        # add comment character:
        if not self._comment : self.line = '!'+self.line
        # reset flag:
        self.setflag( empty=True )
    #enddef
    
    # *
    
    def KeysWithUnresolvedValue( self, marks=('${','}') ) :
    
        """
        If this looks like a 'key : value' line, return the key, 
        or None otherwise. The key should not contain any unresolved
        substitutions, which is tested by presence of the marks,
        and should not start with '#'.
        """
        
        # init result:
        res = {}

        # only if not resolved yet ...
        if not self._resolved :        
            # could this be a 'key : value' line ?
            if ':' in self.line :
                # split, remove leading and end space:
                qkey,qvalue = self.line.split(':',1)
                qkey = qkey.strip()
                # assume it is indeed a key if it does not contain whitespace,
                # no start mark, and does not start wiht '#' :
                if (' ' not in qkey) and (marks[0] not in qkey) and (not qkey.startswith('#')) :
                    ## add to list:
                    #if qkey not in keys_with_unresolved_value : keys_with_unresolved_value.append(qkey)
                    # fill into key/value pair without value:
                    rcv = RcValue( qkey, None, self.line )
                    # store:
                    res = { qkey : rcv }
                #endif
            #endif
        #endif
        
        # ok
        return res
        
    #enddef     # KeysWithUnresolvedValue

    # *
    
    def Extract( self ) :
    
        """
        Extract key/value pair and return {key : RcValue} dictionary.
        """
        
        # no value line or not yet possible to extract key/value ?
        if self._comment or self._special or (not self._resolved) or self._hold :
                
            # empty result:
            res = {}

        else :
        
            # should be a key:value line ..
            if ':' not in self.line :
                logger.error( 'Found a line that is neither empty, comment,' )
                logger.error( 'or a special line, but does not seem a key:value line:' )
                logger.error( self.line )
                logger.error( self.traceback() )
                raise ValueError
            #endif
                    
            # split in key and value; 
            # value might contain ':' too, so at maximum 1 split:
            key,val = self.line.split(':',1)

            # remove comment from value:
            if '!' in val :
                # not if '\!' is in the value ...
                if not '\!' in val : val,comment = val.split('!')
                # replace all slash-comments:
                val = val.replace('\!','!')
            #endif

            # remove spaces:
            key = key.strip()
            val = val.strip()

            # fill into rcvalue:
            rcv = RcValue( key, val, self.copy() )

            # store:
            res = { key : rcv }

        #endif
        
        # ok
        return res

    #enddef

#endclass  # RcLine


# ***


class RcLines( list ) :

    """
    Storage for source lines.
    """
    
    def ReadFile( self, filename, raw=False ) :
    
        """
        Read sourcelines from a file.
        """
        
        # modules:
        import os
        
        # check ...
        if not os.path.exists( filename ) :
            logger.error( 'file not found : %s' % filename )
            raise IOError
        #endif

        # read lines:
        f = open( filename, 'r' )
        lines = f.readlines()
        f.close()

        # loop over lines:
        continuation = False
        for iline in range(len(lines)) :
            # current value, without newline:
            line = lines[iline].rstrip()
            # continuation or new fresh line?
            if continuation :
                # add current line to latest:
                self[len(self)-1].extend( line.rstrip('\\'), iline+1 )
            else :
                # new line:
                newline = RcLine( line.rstrip('\\'), filename, iline+1 )
                # raw mode ? then no need to resolve:
                if raw : newline.setflag( resolved=True )
                # store line and traceback information in new element:
                self.append( newline )
            #endif
            # should be continued with next line ?
            continuation = line.endswith('\\')
        #endfor
        
    #enddef  # ReadFile
    
    # *
    
    def ReadTree( self, SL, upto_mark=None ) :

        """
        Pop lines from RcLines and return a source tree.
        """
        
        # until no lines are left:
        while len(SL) > 0 :
        
            # special line ?
            if SL[0]._special :
                # init flag:
                reached_end = False
                # include line ?
                if SL[0].startswith('#include') :
                    # extract include-block:
                    subtree = RcLinesInclude()
                    subtree.ReadTree( SL )
                # start of if-block ?
                elif SL[0].startswith('#if') :
                    # extract if-block:
                    subtree = RcLinesIf()
                    subtree.ReadTree( SL )
                # part of if-block, just copy:
                elif SL[0].startswith('#else') or \
                     SL[0].startswith('#elif') :
                    # extract from list:
                    subtree = SL.pop(0)
                # error statement ?
                elif SL[0].startswith('#error') :
                    # extract error line:
                    subtree = RcLinesError()
                    subtree.ReadTree( SL )
                # start of for-loop ?
                elif SL[0].startswith('#for') :
                    # extract if-block:
                    subtree = RcLinesFor()
                    subtree.ReadTree( SL )
                # end mark ?
                elif (upto_mark != None) and SL[0].startswith(upto_mark) :
                    # extract from list:
                    subtree = SL.pop(0)
                    # set flag:
                    reached_end = True    
                # other ...
                else :
                    # not yet ...
                    logger.error( 'found unsupported special line:' )
                    logger.error( SL[0] )
                    raise ValueError 
                #endif
                # add to tree:
                self.append( subtree )
                # end was reached ? then leave:
                if reached_end : break

            else :
            
                # pop next line:
                sline = SL.pop(0)
                # add to tree:
                self.append( sline )

            #endif   # special or normal line ?
            
        #endwhile   # lines left

    #enddef
    
    # *
    
    def Show( self, level=0, 
                indent='', autoindent=False, 
                undefined_key=None, marks=('${','}'),
                unresolved=False ) :
    
        """
        Display source lines.
        Indent is displayed at begin of each line.
        If autoindent is set, it is extended with two spaces
        in recursive calls.
        If optional argument 'unresolved' is set to True,
        only unresolved lines are displayed.
        If optional argument 'undefined_key' is set, only
        those lines are displayed that contain an unresolved
        substitution (enclosed by the marks) of this key.
        """
        
        # loop over lines:
        for sline in self :
            # single line ?
            if type(sline)  == RcLine :
                # only if unresolved line ?
                if unresolved and sline.resolved() : continue
                # only if contains the substitution of the given
                # undefined variable ?
                if undefined_key != None :
                    if marks[0]+undefined_key+marks[1] not in sline.getline() : continue
                #endif
                # display:
                logger.log( level, indent+str(sline) )
            else :
                # new indent ?
                if autoindent :
                    newindent = indent+'  '
                else :
                    newindent = indent
                #endif
                # recursive call:
                sline.Show( level=level,
                               indent=newindent, autoindent=autoindent,
                               undefined_key=undefined_key, 
                               unresolved=unresolved )
            #endif
        #endfor
        
    #enddef
    
    # *
    
    def copy( self ) :
        """
        Return copy.
        """
        # copy all elements:
        res = RcLines()
        for i in range(len(self)) :
            res.append( self[i].copy() )
        #endfor
        # ok
        return res
    #enddef
    
    # *
    
    def replace( self, *args ) :
        for i in range(len(self)) : self[i].replace(*args)
    #enddef
    
    # *
    
    def replacevalue( self, *args ) :
        for i in range(len(self)) : self[i].replacevalue( *args )
    #enddef
    
    # *
    
    def setflag( self, **kwargs ) :
    
        """
        Set flags for all sourcelines.
        """
        
        # loop over elements:
        for i in range(len(self)) :
            # set flags:
            self[i].setflag(**kwargs)
        #endif
        
    #enddef
    
    # *
    
    def comment( self ) :
    
        """
        Comment all sourcelines.
        """

        # loop over elements:
        for i in range(len(self)) :
            # comment line(s):
            self[i].comment()
        #endif
        # set flag:
        self._comment = True

    #enddef
    
    # *
    
    def evaluate( self ) :
    
        """
        Try to evaluate sourcelines.
        """
        
        # loop over elements:
        for i in range(len(self)) :
            # evaluate this line or these lines:
            self[i].evaluate()
        #endif
        
    #enddef

    # *
    
    def resolve( self, rcd ) :
    
        """
        Try to resolve sourcelines using dictionary.
        Maintain list with undefined keys.
        """
        
        # init result:
        undefined_keys = []
        
        # loop over elements:
        for i in range(len(self)) :
            # resolve:
            undefkeys = self[i].resolve( rcd )
            # update list:
            for key in undefkeys :
                if key not in undefined_keys : undefined_keys.append(key)
            #endif
        #endif
        
        # ok
        return undefined_keys
        
    #enddef
    
    # *
    
    def resolved( self ) :
        """
        Return True if all elements are resolved.
        """
        # by default true:
        res = True
        # ... unless some are not resolved:
        for sline in self :
            # test for this element:
            res = sline.resolved()
            # not resolved ? then leave:
            if not res : break
        #endfor
        # ok:
        return res
    #enddef  # resolved

    # *
    
    def commented( self ) :
        """
        Return True if all elements are commented.
        """
        # by default true:
        res = True
        # ... unless some are not resolved:
        for sline in self :
            # test for this element:
            res = sline.commented()
            # not commented ? then leave:
            if not res : break
        #endfor
        # ok:
        return res
    #enddef  # commented

    # *
    
    def changed( self ) :
        """
        Return True if some changes were applied after
        the latest call to 'resolved()'.
        """
        # by default false:
        res = False
        # ... unless some are changed:
        for sline in self :
            # test for this element:
            res = sline.changed()
            # was this one changed ? then leave:
            if res : break
        #endfor
        # ok:
        return res
    #enddef  # resolve

    # *
    
    def special( self ) :
        """
        Dummy; always false for collection of lines.
        """
        return False
    #enddef  # resolve

    # *
    
    def Extract( self ) :
    
        """
        Extract key/value pairs and return in a RcValues dictionary.
        """
        
        # empty result:
        res = RcValues()
        # loop over elements:
        for sline in self :
            # get dictionary:
            rcvs = sline.Extract()
            # add values:
            for val in rcvs.values() : res.append( val )
        #endfor
       
        # ok
        return res
        
    #enddef
    
    # *
    
    def KeysWithUnresolvedValue( self, marks=('${','}') ) :
    
        """
        Extract unresolved keys from the lines
        and return in a RcValues dictionary.
        """
        
        # empty result:
        res = RcValues()
        # loop over elements:
        for sl in self :
            # if this line(s) is (contains) key definitions with
            # unresolved values, return a dictionary with the key
            # and an RcValue with None as value:
            rcvs = sl.KeysWithUnresolvedValue(marks=marks)
            # add values:
            for val in rcvs.values() : res.append( val )
        #endfor
       
        # ok
        return res
        
    #enddef
    
    # *
    
    def getlines( self ) :
    
        """
        Return list with lines of rcfile.
        """
        
        # init result:
        res = []
        
        # loop over elements:
        for sl in self :
            # single line or sub-tree ?
            if type(sl) == RcLine :
                # add line:
                res.append( sl.getline() )
            else :
                # recursive call:
                res = res + sl.getlines()
            #endif
        #endfor
        
        # ok:
        return res

    #enddef                
    

#endclass  # RcLines


# ***


class RcLinesInclude( RcLines ) :

    """
    Sourcelines of an '#include ...' line
    """
    
    # *
    
    def ReadTree( self, SL ) :
    
        """
        Extract include line.
        """
        
        # pop first line:
        sline = SL.pop(0)

        # check ...
        if not sline.startswith('#include') :
            logger.error( 'bug: RcLinesInclude.ReadTree should start with "#include"' )
            raise ValueError
        #endif

        # add to (now empty) list:
        self.append( sline )
        
        # internal flags:
        self._evaluated = False

    #enddef
    
    # *
    
    def copy( self ) :
        """
        Return copy.
        """
        # copy all elements:
        res = RcLinesInclude()
        for i in range(len(self)) : res.append( self[i].copy() )
        res._evaluated = self._evaluated
        # ok
        return res
    #enddef
    
    # *

    def evaluate( self ) :
    
        """
        Try to evaluate include line
        """
        
        # leave if all commented ...
        if self.commented() : return

        # not evaluated yet ?
        if (not self._evaluated) :
            # filename is resolved ? then evaluate:
            if self[0].resolved() :
                # remove leading mark, closing colon, and surrounding whitspace:
                includefile = self[0].strip().lstrip('#include').rstrip(':').strip()
                # comment this line:
                self[0].comment()
                # add extra comment:
                sl = self[0].copy()
                sl.setline( '!>>> %s >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>' % includefile )
                self.append( sl )
                # read new sub-tree with inherited routine:
                SL = RcLines()
                SL.ReadFile( includefile )
                # extract sub-tree and add to current list:
                RcLines.ReadTree( self, SL )
                # add extra comment:
                sl = self[0].copy()
                sl.setline( '!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<' )
                self.append( sl )
                # reset flag:
                self._evaluated = True
            #endif
        #endif   # evaluated

        # loop over elements:
        for i in range(len(self)) :
            # evaluate this line or these lines:
            self[i].evaluate()
        #endif
        
    #enddef    # evaluate
    
#endclass   # RcLinesInclude


# ***


class RcLinesIf( RcLines ) :

    """
    Sourcelines of an '#if ... #endif' block.
    """
    
    # *
    
    def ReadTree( self, SL ) :
    
        """
        Extract if-endif block.
        """
        
        # pop first line:
        sline = SL.pop(0)

        # check ...
        if not sline.startswith('#if') :
            logger.error( 'bug: RcLinesIf.ReadTree should start with "#if"' )
            raise ValueError
        #endif

        # extract new sub-tree with inherited routine:
        RcLines.ReadTree( self, SL, upto_mark='#endif' )
        # prepend start line:
        self.insert(0,sline)
        
        # check ...
        nline = len(self)
        if not self[nline-1].startswith('#endif') :
            logger.error( '#if without #endif at %s' % self[0].traceback() )
            raise ValueError
        #endif
        
        # hold all lines except first if-statement
        for i in range(1,nline) : self[i].setflag( hold=True )
        
        # internal data:
        self._iline = 0   # condition line: leading #if or #elif

    #enddef
    
    # *
    
    def copy( self ) :
        """
        Return copy.
        """
        # copy all elements:
        res = RcLinesIf()
        for i in range(len(self)) : res.append( self[i].copy() )
        res._iline = self._iline
        # ok
        return res
    #enddef
    
    # *
    
    def setflag( self, **kwargs ) :
    
        """
        Set flags for all lines or only up to condition line.
        """

        # loop over elements:
        for i in range(len(self)) :
            # condition to be evaluated ? then only until that line:
            if (self._iline != None) and (i > self._iline) : break
            # set flags:
            self[i].setflag(**kwargs)
        #endif
        
    #enddef

    # *

    def evaluate( self ) :
    
        """
        Try to evaluate if-then-else condition.
        """
        
        # leave if all commented ...
        if self.commented() : return

        ## debug ...
        #logger.info( '  evaluate if-block ; iline=%s' % str(self._iline) )
        
        # line to be tested ?
        if self._iline != None :
            # condition line is resolved ? then evaluate:
            if self[self._iline].resolved() :
                # remove leading mark, closing colon, and surrounding whitspace:
                condition = self[self._iline].strip().lstrip('#if').lstrip('#elif').rstrip(':').strip()
                # evaluate the condition:
                flag = eval( condition )
                # flag to denote if the current line is evaluated:
                is_evaluated = flag
                # comment #if line:
                self[self._iline].comment()
                # save index:
                iline = self._iline
                # reset line index:
                self._iline = None
                # loop over remaining lines :
                search_first = True
                for i in range(iline+1,len(self)) :
                    # no special line found yet ?
                    if search_first :
                        # special line ?
                        if self[i].special() :
                            # which special ?
                            if self[i].startswith('#elif') :
                                # swap flag to (un)comment lines inbetween:
                                flag = not flag
                                # is this condition to be tested now ?
                                if flag :
                                    # store line number:
                                    self._iline = i
                                    # unhold this line:
                                    self[i].setflag( hold=False, changed=True )
                                    # remainder should be evaluated in a next pass:
                                    is_evaluated = False
                                #endif
                            elif self[i].startswith('#else') :
                                # this seems to be a #if..#else..#endif block;
                                # the #else line where we are now could be commented, 
                                # as well as the closing #endif line:
                                self[i].comment()
                                self[len(self)-1].comment()
                                # swap flag to (un)comment lines in between:
                                flag = not flag
                                # reset flag:
                                is_evaluated = flag
                                # next:
                                continue
                            elif self[i].startswith('#endif') :
                                # this seems to be a simple #if..#endif
                                # block which is now evaluated, so this
                                # is the final #endif line which could be
                                # commented now:
                                self[i].comment()
                                # jump out:
                                break
                            else :
                                # something wrong ...
                                logger.error( 'unsupported special line: %s' % self[i] )
                                raise ValueError
                            #endif
                            # reset flag:
                            search_first = False
                        #endif
                    #endif  # searching first special
                    # is this a True or a False part ?
                    if flag :
                        # True part could now be resolved:
                        if is_evaluated : self[i].setflag( hold=False, changed=True )
                    else :
                        # False part, so comment this line:
                        self[i].comment()
                    #endif
                #endfor   # lines
            #endif   # condition line resolved ?
        #endif   # condition to be evaluated
                
        # loop over elements:
        for i in range(len(self)) :
            # evaluate this line or these lines:
            self[i].evaluate()
        #endif
        
    #enddef    # evaluate
    
#endclass   # RcLinesIf


# ***


class RcLinesError( RcLines ) :

    """
    Sourceline of an '#error ...' statement.
    """
    
    # *
    
    def ReadTree( self, SL ) :
    
        """
        Extract error line.
        """
        
        # pop first line:
        sline = SL.pop(0)

        # check ...
        if not sline.startswith('#error') :
            logger.error( 'bug: RcLinesInclude.ReadTree should start with "#error"' )
            raise ValueError
        #endif

        # add to (now empty) list:
        self.append( sline )
        
    #enddef
    
    # *
    
    def copy( self ) :
        """
        Return copy.
        """
        # copy all elements:
        res = RcLinesError()
        for i in range(len(self)) : res.append( self[i].copy() )
        # ok
        return res
    #enddef
    
    # *

    def evaluate( self ) :
    
        """
        Try to evaluate include line
        """
        
        # leave if all commented ...
        if self.commented() : return

        # message is not empty and resolved ? then evaluate:
        if self[0].resolved() :
            # remove leading mark, closing colon, and surrounding whitspace:
            msg = self[0].strip().lstrip('#error').rstrip(':').strip()
            # display:
            logger.error( msg )
            # add info:
            logger.error( 'error message in %s' % self[0].traceback() )
            # stop:
            raise Exception
        #endif
        
    #enddef    # evaluate
    
#endclass   # RcLinesError


# ***


class RcLinesFor( RcLines ) :

    """
    Sourcelines of an '#for ... #endfor' block.
    """
    
    # *
    
    def ReadTree( self, SL ) :
    
        """
        Extract for-endif block.
        """
        
        # pop first line:
        sline = SL.pop(0)

        # check ...
        if not sline.startswith('#for') :
            logger.error( 'bug: RcLinesFor.ReadTree should start with "#for"' )
            raise ValueError
        #endif

        # extract new sub-tree with inherited routine:
        RcLines.ReadTree( self, SL, upto_mark='#endfor' )
        # prepend start line:
        self.insert(0,sline)
        
        # check ...
        nline = len(self)
        if not self[nline-1].startswith('#endfor') :
            logger.error( '#for without #endfor at %s' % self[0].traceback() )
            raise ValueError
        #endif

        # hold all lines except first for-statement:
        for i in range(1,nline) : self[i].setflag( hold=True )
        
        # internal flags:
        self._evaluated = False

    #enddef

    # *
    
    def copy( self ) :
        """
        Return copy.
        """
        # copy all elements:
        res = RcLinesFor()
        for i in range(len(self)) : res.append( self[i].copy() )
        res._evaluated = self._evaluated
        # ok
        return res
    #enddef
    
    # *

    def evaluate( self ) :
    
        """
        Try to evaluate if-then-else condition.
        """
        
        # leave if all commented ...
        if self.commented() : return

        # for-line not evaluated yet ?
        if not self._evaluated :
            # for-line is resolved ? then evaluate:
            if self[0].resolved() :
                # remove leading mark, closing colon, and surrounding whitspace:
                expression = self[0].strip().lstrip('#for').rstrip(':').strip()
                # split at whitspace:
                elements = expression.split()
                # check for at least:  ["var","in"]
                if (len(elements) < 2) or (elements[1] != 'in') :
                    logger.error( 'for-loop expression should be of form "for var in"' )
                    logger.error( str(self[0]) )
                    raise ValueError
                #endif
                # extract variable name and values:
                varname = elements[0]
                values = elements[2:]
                # comment the #for and the #endfor lines:
                self[0].comment()
                self[len(self)-1].comment()
                # extract lines:
                slines = []
                for i in range(len(self)-2) : slines.append( self.pop(1) )
                # loop over elements:
                for value in values :
                    # insert lines:
                    for sline in slines :
                        # new source line:
                        sl = sline.copy()
                        # replace variable with value:
                        sl.replace(varname,value)
                        # unhold:
                        sl.setflag( hold=False )
                        # insert new source line:
                        self.insert( len(self)-1, sl )
                    #endfor
                #endfor
                # reset flag:
                self._evaluated = True
            #endif   # condition line resolved ?
        #endif   # condition to be evaluated
                
        # loop over elements:
        for i in range(len(self)) :
            # evaluate this line or these lines:
            self[i].evaluate()
        #endif
        
    #enddef    # evaluate
    
#endclass   # RcLinesFor


#-------------------------------------------------
# rcvalue(s)
#-------------------------------------------------


class RcValue( object ) :

    """
    Storage for key/value pair and traceback info.
    """
    
    def __init__( self, key, value, sourceline ) :
        """
        Store rcvalue data.
        """
        # store:
        self.key = key
        self.value = value
        self.sourceline = sourceline
    #enddef
    
    # *
    
    def getkey( self ) :
        """
        Return rc key.
        """
        return self.key
    #enddef
    
    # *
    
    def getvalue( self ) :
        """
        Return rc value.
        """
        return self.value
    #enddef
    
    
    # *
    
    def getsourceline( self ) :
        """
        Return rc sourceline.
        """
        return self.sourceline
    #enddef
    
    # *
    
    def replacevalue( self, value ) :
        """
        Replace value, and also change source line.
        """
        # store new value:
        self.value = value
        # replace source line:
        self.sourceline.setline( '%s : %s' % (self.key,value) )
    #enddef
    
#endclass


# ***


class RcValues( dict ) :

    """
    Storage for key/value pairs and traceback info.
    """
    
    def append( self, rcv ) :
    
        """
        Add RcValue to dictionary.
        """
        
        # already present ?
        if rcv.getkey() in self.keys() :
            logger.error( 'key "%s" defined at least twice :' % rcv.getkey() )
            logger.error( self[rcv.getkey()].getsourceline() )
            logger.error( rcv.getsourceline() )
            raise ValueError
        #endif
        
        # add:
        self[rcv.key] = rcv
        
    #enddef

#endclass



#-------------------------------------------------
# routines (for backward compatibility)
#-------------------------------------------------


def read( rcfilename, silent=False ) :

    """ 
    This method reads an rc-file by making an instance of the RcFile class, 
    and then returns the dictionary of values only. 
    This makes it backwards compatible with older implementations of the rc.py module
    """

    # read file into RcFile object:
    rcf = RcFile( rcfilename, silent=silent )
    # fill result:
    res = {}
    for key in rcf.keys() : res[key] = rcf.get(key)
    
    # ok
    return res

#enddef


# ***


def write( filename, rcdict ) :

    """
    This method writes an rc-file dictionary. 
    This makes it backwards compatible with older implementations of the rc.py module
    """

    # open file for writing:
    f = open(filename,'w')

    # loop over key/value pairs:
    for k,v in rcdict.items():
        # add line; at least the specified number of characters 
        # is used for the key:
        f.write( '%-20s:%s\n' % (k,v) )
    #endfor

    # close file:
    f.close()

#enddef


#-------------------------------------------------
# end
#-------------------------------------------------
