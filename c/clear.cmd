@echo Removing temporary files

@del *.resharper.user
@rem del *.resharper
@del *.cache
rem @del /A:H *.suo
@rmdir /S/Q _ReSharper.ExactFFT
@rmdir /S/Q TestResults
@del *.xml
@del *.user
@del *.cgl
@del *.obj
@del *.tds
@del *.sdf
@del *.~*
@rmdir /S/Q bin
@rmdir /S/Q obj
@rmdir /S/Q Debug
@rmdir /S/Q Release