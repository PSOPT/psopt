Const ForReading = 1
Const ForWriting = 2

Set objFSO = CreateObject("Scripting.FileSystemObject")
Set objTextFile = objFSO.OpenTextFile(".\lusol\csrc\myblas.h", ForReading)

Do Until objTextFile.AtEndOfStream
    strNextLine = objTextFile.Readline

    intLineFinder = InStr(strNextLine, "windows.h")
    If intLineFinder <> 0 Then
        strNextLine = "// #include <windows.h>"
    End If

    strNewFile = strNewFile & strNextLine & vbCrLf
Loop

objTextFile.Close

Set objTextFile = objFSO.OpenTextFile(".\lusol\csrc\myblas.h", ForWriting)

objTextFile.WriteLine strNewFile
objTextFile.Close

