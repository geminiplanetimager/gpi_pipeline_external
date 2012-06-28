;+
; NAME:
;       ErrorLogger__Define
;
; PURPOSE:
;
;       The purpose of this program is to log program errors or text messages during
;       program execution as an aid it debugging such a program at a later date. The
;       ErrorLogger program is written as an object so that it will persist in the IDL
;       session until it is destroyed.
;
; AUTHOR:
;
;       FANNING SOFTWARE CONSULTING
;       David Fanning, Ph.D.
;       1645 Sheely Drive
;       Fort Collins, CO 80526 USA
;       Phone: 970-221-0438
;       E-mail: davidf@dfanning.com
;       Coyote's Guide to IDL Programming: http://www.dfanning.com
;
; CATEGORY:
;
;       Utilities
;
; CALLING SEQUENCE:
;
;       errorLogger = Obj_New("ErrorLogger")
;
; ARGUMENTS:
;
;       filename:    The name of the error log file. If not provided, a default name
;                    will be created, based on the current system time. (Optional)
;
; KEYWORDS:
;
;       ALERT:       The default behavior of the error logger is simply to write text to a file.
;                    But if the ALERT keyword is set, the program will alert the user via a
;                    message dialog that an error has occurred when using the AddError method. (Input)
;
;      DELETE_ON_DESTROY: If this keyword is set, the error log file will be deleted when the
;                    ErrorLogger object is destroyed, but only if the ErrorLogger object is not
;                    in an error state at that time.
;
;       NOTRACEBACK: Set this keyword to suppress traceback information in the error log output
;                    and in any alerts issued by the program. (Input)
;
;       TIMESTAMP:   Set this keyword if you wish a time stamp to be appended to the provided
;                    filename. Otherwise, the filename is used as defined. Default filenames
;                    always have a timestamp appended to the file name. (Input)
;
; METHODS:
;
;        AddError:   Adds an error text string or array to the error log file. By default,
;                    it will add the HELP, LAST_MESSAGE=1, /TRACEBACE traceback 
;                    information to the file. (Procedure)
;
;        AddText:    Adds a text string or array to the error log file. (Procedure)
;
;        ClearLog:   Erases all the text currently in the error log file. (Procedure)
;
;        CloseFile:  Closes the currently open error log file. (Procedure)
;
;        GetProperty: Gets properties of the object. (Procedure)
;
;        LastMessage: Returns the last message text written into the error log file. (Function)
;
;        OpenFile:   Opens the error log file for writing. (Function)
;
;        PrintLastMessage: Writes the last message text written into the error log file to 
;                    standard output. (Procedure)
;
;        Status:     Returns the current status of the error logger. (0 - waiting for input, 
;                    1 - normal operation, 2 - error operation.) (Function)
;
;        SetProperty: Sets properties of the object. (Procedure)
;
;        SetStatus:  Sets the current status of the error logger. Normally not used by the
;                    user, but used internally. (Procedure)
;
; MODIFICATION HISTORY:
;
;       Written by David W. Fanning, November 2009.
;       Modified and expanded the way errors are written into the log file and displayed.
;          Also made it possible to automatically delete the log file when the object is
;          destroyed, if the error logger is not in an error state at the time. Added
;          DELETE_ON_DESTROY and NOTRACEBACK keywords to the INIT and SetProperty
;          methods. 28 Jan 2010. DWF.
;        Modified default filenames so that I am now guaranteed to get unique file names 
;           by using Timestamp program from the Coyote Library. 8 Feb 2010. DWF
;-
;
;******************************************************************************************;
;  Copyright (c) 2009, by Fanning Software Consulting, Inc.                                ;
;  All rights reserved.                                                                    ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;
;+
; NAME:
;       ErrorLogger::AddError
;
; PURPOSE:
;
;       Adds error text to the error log file and sets the error log status to 2 (error
;       condition). If the error logger alert flag is set to 1, the method will alert
;       the user to the error with a pop-up message dialog as well as writing the output
;       to standard output.
;
; CALLING SEQUENCE:
;
;       errorLogger -> AddError, theText
;
; ARGUMENTS:
;
;       theText :    The error message text you wish to add to the file. If not provided,
;                    the text of the last error message (Help, /LAST_MESSAGE) is used and
;                    written to the file.
;
; KEYWORDS:
;
;       None.
;-
;******************************************************************************************;
PRO ErrorLogger::AddError, theText

   ; The text will equal the ERROR_STATE messsage string if not available.
   IF N_Elements(theText) EQ 0 THEN theText = !Error_State.Msg
   
   ; Get the call stack and the calling routine's name.
   callStack = Scope_Traceback()
   callingRoutine = (StrSplit(StrCompress(callStack[1])," ", /Extract))[0]
   
   ; Are widgets supported?
   IF !D.Name EQ 'PS' THEN BEGIN
      widgetsSupported = 1
   ENDIF ELSE BEGIN
      widgetsSupported = ((!D.Flags AND 65536L) NE 0)
   ENDELSE

   ; No dialogs, unless alert is set.
   IF ~Keyword_Set(self.alert) THEN widgetsSupported = 0
   
   ; It is not enough to know if widgets are supported. In CRON jobs, widgets are
   ; supported, but there is no X connection and pop-up dialogs are not allowed.
   ; Here is a quick test to see if we can connect to a windowing system. If not,
   ; then we are going to assume widgets are not supported.
   Catch, theError
   IF theError NE 0 THEN BEGIN
      Catch, /CANCEL
      widgetsSupported = 0
      GOTO, testWidgetSupport
   ENDIF
   theWindow = !D.Window
   IF (!D.Flags AND 256) NE 0 THEN Window, /FREE, XSIZE=5, YSIZE=5, /PIXMAP
   Catch, /CANCEL
   
   testWidgetSupport: ; Come here if you choke on creating a window.
   IF !D.Window NE theWindow THEN BEGIN
      WDelete, !D.Window
      IF theWindow GE 0 THEN WSet, theWindow
   ENDIF
   
   IF widgetsSupported && self.alert THEN BEGIN
   
      ; If this is an error produced with the MESSAGE command, it is a trapped
      ; error and will have the name "IDL_M_USER_ERR".
      IF !ERROR_STATE.NAME EQ "IDL_M_USER_ERR" THEN BEGIN
   
         IF N_Elements(title) EQ 0 THEN title = 'Trapped Error'
   
            ; If the message has the name of the calling routine in it,
            ; it should be stripped out. Can you find a colon in the string?
   
         ; Is the calling routine an object method? If so, special processing
         ; is required. Object methods will have two colons together.
         doublecolon = StrPos(theText, "::")
         IF doublecolon NE -1 THEN BEGIN
   
            prefix = StrMid(theText, 0, doublecolon+2)
            submessage = StrMid(theText, doublecolon+2)
            colon = StrPos(submessage, ":")
            IF colon NE -1 THEN BEGIN
   
               ; Extract the text up to the colon. Is this the same as
               ; the callingRoutine? If so, strip it.
               IF StrMid(theText, 0, colon+StrLen(prefix)) EQ callingRoutine THEN $
                  theText = StrMid(theText, colon+1+StrLen(prefix))
            ENDIF
         ENDIF ELSE BEGIN
   
            colon = StrPos(theText, ":")
            IF colon NE -1 THEN BEGIN
   
               ; Extract the text up to the colon. Is this the same as
               ; the callingRoutine? If so, strip it.
               IF StrMid(theText, 0, colon) EQ callingRoutine THEN $
                  theText = StrMid(theText, colon+1)
            ENDIF
   
         ENDELSE
   
         ; Add the calling routine's name.
         void = Dialog_Message(StrUpCase(callingRoutine) + ": " + theText, Title=title)
   
      ENDIF ELSE BEGIN
   
         ; Otherwise, this is an IDL system error.
         IF N_Elements(title) EQ 0 THEN title = 'System Error'
   
         IF StrUpCase(callingRoutine) EQ "$MAIN$" THEN BEGIN
            void = Dialog_Message(theText, Title=title) 
         ENDIF ELSE BEGIN
             void = Dialog_Message(StrUpCase(callingRoutine) + "--> " + theText, Title=title)
         ENDELSE
      ENDELSE
   ENDIF ELSE BEGIN
         Message, theText, /Continue, /NoPrint, /NoName, /NoPrefix
         self -> AddText, '%' + callingRoutine + ': ' + theText 
   ENDELSE
   
   ; Provide traceback information if requested.
   Help, /Last_Message, Output=traceback
   traceback = ['Traceback Report from ' + StrUpCase(callingRoutine) + ':', $
                 '', "     " + traceback]
   IF self.alert THEN BEGIN
       IF Keyword_Set(self.traceback) THEN BEGIN
          FOR j=0,N_Elements(traceback)-1 DO Print, "     " + traceback[j]
       ENDIF
   ENDIF 
   IF Keyword_Set(self.traceback) THEN self -> AddText, traceback

    ; Set the status to error condition.
    self -> SetStatus, 2
    
END


;******************************************************************************************;
;+
; NAME:
;       ErrorLogger::AddText
;
; PURPOSE:
;
;       Adds text to the error log file and sets the error log status to 1 (normal
;       condition). 
;
; CALLING SEQUENCE:
;
;       errorLogger -> AddText, theText
;
; ARGUMENTS:
;
;       theText :    The message text you wish to add to the file. 
;
; KEYWORDS:
;
;       None.
;-
;******************************************************************************************;
PRO ErrorLogger::AddText, theText
    
    Compile_Opt idl2
    
    ; Error handling.
    Catch, theError
    IF theError NE 0 THEN BEGIN
        Catch, /CANCEL
        void = Error_Message()
        RETURN
    ENDIF
    
    ; Have to have text to do anything.
    IF N_Elements(theText) EQ 0 THEN RETURN
    
    ; Make sure these are strings we are writing.
    thisType = Size(theText, /TNAME)
    IF thisType NE 'STRING' THEN Message, 'Only strings can be written into the error log file.'
    
    ; Write the text to the file.
    numLines = N_Elements(theText)
    PrintF, self.lun, ""
    FOR j=0L, N_Elements(theText) -1 DO BEGIN
        PrintF, self.lun, theText[j]
    ENDFOR
    
    ; Update the error logger status to normal.
    self.status = 1
    
    ; Save the last message for later recall.
    *self.lastMessage = theText
    
    
END 


;******************************************************************************************;
;+
; NAME:
;       ErrorLogger::ClearLog
;
; PURPOSE:
;
;       Clears the error log file of text.
;
; CALLING SEQUENCE:
;
;       errorLogger -> ClearLog
;
; ARGUMENTS:
;
;       None.
;
; KEYWORDS:
;
;       None.
;-
;******************************************************************************************;
PRO ErrorLogger::ClearLog

    ; Close the current error log file and delete it.
    self -> CloseFile
    File_Delete, self.filename, /ALLOW_NONEXISTENT
    
    ; Open a new error log file with the same name.
    self -> OpenFile, self.filename
    
    ; Set the error logger status to waiting.
    self.status = 0
    
END 


;******************************************************************************************;
;+
; NAME:
;       ErrorLogger::CloseFile
;
; PURPOSE:
;
;       Closes the currently open error log file.
;
; CALLING SEQUENCE:
;
;       errorLogger -> CloseFile
;
; ARGUMENTS:
;
;       None.
;
; KEYWORDS:
;
;       None.
;-
;******************************************************************************************;
PRO ErrorLogger::CloseFile
    IF self.lun GE 100 THEN Free_Lun, self.lun ELSE IF self.lun GT 0 THEN Close, self.lun
END 


;******************************************************************************************;
;+
; NAME:
;       ErrorLogger::GetFileName
;
; PURPOSE:
;
;       Returns the file name of the error log file.
;
; CALLING SEQUENCE:
;
;       filename =  errorLogger -> GetFileName()
;
; RETURN VALUE:
;
;       filename:     The name of the error log file.
; 
; ARGUMENTS:
;
;       None.
;
; KEYWORDS:
;
;       None.
;-
;******************************************************************************************;
FUNCTION ErrorLogger::GetFileName
    RETURN, self.filename
END 



;******************************************************************************************;
;+
; NAME:
;       ErrorLogger::Status
;
; PURPOSE:
;
;       Returns the current status of the error logger.
;
; CALLING SEQUENCE:
;
;       status = errorLogger -> Status()
;
; RETURN VALUE:
;
;       status:     The error log status: 
;                      0 - waiting for input
;                      1 - normal operation
;                      2 - error operation
;
; ARGUMENTS:
;
;       None.
;
; KEYWORDS:
;
;       None.
;-
;******************************************************************************************;
FUNCTION ErrorLogger::Status
    RETURN, self.status
END 


;******************************************************************************************;
;+
; NAME:
;       ErrorLogger::SetStatus
;
; PURPOSE:
;
;       Sets the current status of the error logger.
;
; CALLING SEQUENCE:
;
;       errorLogger -> SetStatus, status
;
; ARGUMENTS:
;
;       status:     The error log status: 
;                      0 - waiting for input
;                      1 - normal operation
;                      2 - error operation
;
; KEYWORDS:
;
;       None.
;-
;******************************************************************************************;
PRO ErrorLogger::SetStatus, status
    IF N_Elements(status) NE 0 THEN self.status = 0 > status < 2
END 



;******************************************************************************************;
;+
; NAME:
;       ErrorLogger::OpenFile
;
; PURPOSE:
;
;       Opens the error log file.
;
; CALLING SEQUENCE:
;
;       errorLogger -> OpenFile, filename
;
; ARGUMENTS:
;
;       filename:     The name of the error log file.
; 
; KEYWORDS:
;
;       DELETE_CURRENT_FILE:  If this keyword is set, the current error log file is closed
;                      and deleted before the new file is opened for writing.
;-
;******************************************************************************************;
FUNCTION ErrorLogger::OpenFile, newLogFilename, DELETE_CURRENT_FILE=delete_current_file
    
    Compile_Opt idl2
    
    ; Error handling
    Catch, theError
    IF theError NE 0 THEN BEGIN
        Catch, /CANCEL
        void = Error_Message()
        IF N_Elements(lun) NE 0 THEN BEGIN
            Free_Lun, lun
            File_Delete, newLogFilename, /ALLOW_NONEXISTENT
        ENDIF
        RETURN, 0
    ENDIF
    
    ; Can we write into the specified directory?
    basename = FSC_Base_Filename(newLogFilename, EXTENSION=ext, DIRECTORY=dir)
    IF File_Test(dir, /DIRECTORY) EQ 0 THEN Message, 'Specified directory (' + dir + ') does not exist.' 
    
    ; Close the current file (if any) before opening a new one.
    self -> CloseFile
    
    ; Need to delete the current file?
    IF Keyword_Set(delete_current_file) THEN File_Delete, self.filename, /ALLOW_NONEXISTENT

    ; Open the file for writing.
    OpenW, lun, newLogFilename, /GET_LUN
    self.lun = lun
    
    ; Write a header into the file.
    PrintF, self.lun, 'Error log file created: ' + Systime()
    PrintF, self.lun, ""
    
    ; Store the filename
    self.filename = newLogFilename
    
    RETURN, 1
    
END


;******************************************************************************************;
;+
; NAME:
;       ErrorLogger::LastMessage
;
; PURPOSE:
;
;       Returns the last text message written to the error logger.
;
; CALLING SEQUENCE:
;
;       message = errorLogger -> LastMessage()
;
; RETURN VALUE:
;
;       message:     The last text written to the error log file.
;
; ARGUMENTS:
;
;       None.
;
; KEYWORDS:
;
;       None.
;-
;******************************************************************************************;
FUNCTION ErrorLogger::LastMessage
    
    ; Returns the last message added to the file.
    IF N_Elements(*self.lastMessage) NE 0 THEN RETURN, *self.lastMessage ELSE RETURN, ""
    
END 


;******************************************************************************************;
;+
; NAME:
;       ErrorLogger::PrintLastMessage
;
; PURPOSE:
;
;       Prints the last text message written to the error logger to standard output.
;
; CALLING SEQUENCE:
;
;       errorLogger -> PrintLastMessage
;
; ARGUMENTS:
;
;       None.
;
; KEYWORDS:
;
;       None.
;-
;******************************************************************************************;
PRO ErrorLogger::PrintLastMessage
    
    ; Prints the last message in the error logger.
    lastMessage = self -> LastMessage()
    FOR j=0, N_Elements(lastMessage)-1 DO Print, lastMessage[j]
    
END 



;******************************************************************************************;
;+
; NAME:
;       ErrorLogger::GetProperty
;
; PURPOSE:
;
;       Allows the user to get properties from the object via keywords.
;
; CALLING SEQUENCE:
;
;       errorLogger -> GetProperty, ...
;
; ARGUMENTS:
;
;       None.
;
; KEYWORDS:
;
;       ALERT:          The alert flag in the object. (Output)
;
;       DELETE_ON_DESTROY:  The delete on destroy flag in the object. (Output)
;       
;       FILENAME:       The name of the error log file. (Output)
;
;       LAST_MESSAGE:   The last message written into the error log file. (Output)
;
;       LUN:            The logical unit number of the open error log file. (Output)
;
;       STATUS:         The current error log status. (Output)
;
;       NOTRACEBACK:    The notraceback flag in the object. (Output)
;
;-
;******************************************************************************************;
PRO ErrorLogger::GetProperty, $
    ALERT=alert, $
    DELETE_ON_DESTROY=delete_on_destroy, $
    FILENAME=filename, $
    LAST_MESSAGE=last_message, $
    LUN=lun, $
    NOTRACEBACK=notraceback, $
    STATUS=status
    
    alert = self.alert
    delete_on_destroy = self.delete_on_destroy
    filename = self.filename
    last_message = self -> LastMessage()
    lun = self.lun
    notraceback = ~self.traceback
    status = self.status

END 



;******************************************************************************************;
;+
; NAME:
;       ErrorLogger::SetProperty
;
; PURPOSE:
;
;       Allows the user to set properties of the object via keywords.
;
; CALLING SEQUENCE:
;
;       errorLogger -> SetProperty, ...
;
; ARGUMENTS:
;
;       None.
;
; KEYWORDS:
;
;       ALERT:          The alert flag in the object. (Input)
;       
;       DELETE_ON_DESTROY:  The delete on destroy flag in the object. (Input)
;       
;       NOTRACEBACK:    The notraceback flag in the object. (Input)
;
;       STATUS:         The current error log status. (Input)
;-
;******************************************************************************************;
PRO ErrorLogger::SetProperty, $
    ALERT=alert, $
    DELETE_ON_DESTROY=delete_on_destroy, $
    NOTRACEBACK=notraceback, $
    STATUS=status
    
    IF N_Elements(alert) NE 0 THEN self.alert = Keyword_Set(alert)
    IF N_Elements(status) NE 0 THEN self.status = 0 > status < 2
    IF N_Elements(delete_on_destroy) NE 0 THEN $
        self.delete_on_destroy = Keyword_Set(delete_on_destroy)
    IF N_Elements(notraceback) NE 0 THEN self.tracebace = 1 - Keyword_Set(notraceback)
    IF N_Elements(status) NE 0 THEN self.status = 0 > status < 2

END 



;******************************************************************************************;
;+
; NAME:
;       ErrorLogger::CLEANUP
;
; PURPOSE:
;
;       Cleans up the object.
;
; CALLING SEQUENCE:
;
;       Called automatically when the object is destroyed.
;
; ARGUMENTS:
;
;       None.
;
; KEYWORDS:
;
;       None.
;-
;******************************************************************************************;
PRO ErrorLogger::CLEANUP

    ; Be sure the file is close. Otherwise, it can't be deleted.
    IF self.lun GE 100 THEN Free_Lun, self.lun ELSE IF self.lun GT 0 THEN Close, self.lun
    
    ; If the file is not in an error state, and the delete_on_destroy flag
    ; is set, delete the error log file.
    IF self.delete_on_destroy THEN BEGIN
        IF self.status NE 2 THEN File_Delete, self.filename, /ALLOW_NONEXISTENT
    ENDIF
    
    ; Free the last message pointer.
    Ptr_Free, self.lastMessage
    
END 


;******************************************************************************************;
;+
; NAME:
;       ErrorLogger::INIT
;
; PURPOSE:
;
;       The initialization method for the object.
;
; CALLING SEQUENCE:
;
;       check = Obj_New('ErrorLogger', filename)
;
; ARGUMENTS:
;
;       filename:    The name of the error log file. If not provided, a default name
;                    will be created based on the current system time. (Optional)
;
; KEYWORDS:
;
;       ALERT:       The default behavior of the error logger is simply to write text to a file.
;                    But if the ALERT keyword is set, the program will alert the user via a
;                    message dialog that an error has occurred when using the AddError method. (Input)
;                    
;      DELETE_ON_DESTROY: If this keyword is set, the error log file will be deleted when the
;                    ErrorLogger object is destroyed, but only if the ErrorLogger object is not
;                    in an error state at that time.
;
;       NOTRACEBACK: Set this keyword to suppress traceback information in the error log output
;                    and in any alerts issued by the program. (Input)
;
;       TIMESTAMP:   Set this keyword if you wish a time stamp to be appended to the provided
;                    filename. Otherwise, the filename is used as defined. Default filenames
;                    always have a timestamp appended to the file name. (Input)
;                    
;-
;******************************************************************************************;
FUNCTION ErrorLogger::INIT, filename, $
    ALERT=alert, $
    DELETE_ON_DESTROY=delete_on_destroy, $
    NOTRACEBACK=notraceback, $
    TIMESTAMP=timestamp
    
    COMPILE_OPT idl2

    Catch, theError
    IF theError NE 0 THEN BEGIN
        Catch, /CANCEL
        void = Error_Message()
        IF N_Elements(lun) NE 0 THEN BEGIN
            Free_Lun, lun
            File_Delete, logFilename, /ALLOW_NONEXISTENT
        ENDIF
        RETURN, 0
    ENDIF

    ; Does the filename exist?
    IF N_Elements(filename) EQ 0 THEN BEGIN
       CD, CURRENT=currentDir
       logFilename = FilePath(ROOT_DIR=currentDir, 'logger' + $
            Timestamp(RANDOM_DIGITS=6, /VALID) + '.log')
       filename = logFilename
       timestamp = 0
    ENDIF

    ; Is this a fully qualified filename?
    baseName = File_BaseName(filename)
    IF baseName EQ filename THEN BEGIN
        CD, CURRENT=currentDir
        logFilename = FilePath(ROOT_DIR=currentDir, filename)
    ENDIF ELSE logFilename = filename
    
    ; Does the name need a time stamp?
    IF Keyword_Set(timestamp) THEN BEGIN
       basename = FSC_Base_Filename(logFilename, EXTENSION=ext, DIRECTORY=dir)
       time = Systime(1)
       randomdigits =  StrMid(StrTrim(time - Long(time),2), 2)
       logFilename = Filepath(ROOT_DIR=dir, basename +  Timestamp(RANDOM_DIGITS=6, /VALID))
       IF ext NE "" THEN logFilename = logFilename + '.' + ext
    END

    ; Open the file.
    success = self -> OpenFile(logFilename)
    IF ~success THEN RETURN, 0
    
    ; Store the filename.
    self.filename = logFilename
    
    ; Initialize the last message pointer.
    self.lastMessage = Ptr_New(/ALLOCATE_HEAP)
    
    ; Need alerts?
    self.alert = Keyword_Set(alert)
    
    ; Need traceback?
    self.traceback = 1 - Keyword_Set(notraceback)
    
    ; Delete file on destroy?
    self.delete_on_destroy = Keyword_Set(delete_on_destroy)
    
    ; Successful completion.
    RETURN, 1
    
END 


PRO ErrorLogger__Define, class

   class = { ERRORLOGGER, $
             filename: "", $            ; The error log filename.
             lun: 0L, $                 ; The file logical unit number.
             alert: 0L, $               ; A flag, if set, will give user alerts on errors.
             traceback: 0L, $           ; If set, will include traceback information into the log file.
             lastMessage: Ptr_New(), $  ; The last message written into the file.
             delete_on_destroy: 0L, $   ; A flag causing log file to be deleted when object is destroyed.
             status: 0L }               ; The current status of the error logger. 0-waiting, 1-normal, 2-error.
             
END