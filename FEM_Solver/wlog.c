/* ------------------------------------------------------------------------------------------ */
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include "constants.h"
#include "wlog.h"
/* ------------------------------------------------------------------------------------------ */
char _wlog_file[FILENAME_LENGTH];
/* ------------------------------------------------------------------------------------------ */
void setlogsilent(const char fname[]);
/* ------------------------------------------------------------------------------------------ */
void writelog(const char format[],...)
{
    FILE    *log;
    va_list  argptr;

    if (strlen(_wlog_file) == 0)
    {
        setlogsilent("writelog.log");
    }
    
    log = fopen(_wlog_file, "a+"); 
    
    if (log != NULL)
    {
        va_start(argptr, format);
        vfprintf(log, format, argptr);  
        va_end(argptr);

        va_start(argptr, format);
        vfprintf(stdout, format, argptr);
        va_end(argptr);
        
        fclose(log);
    }
}
/* ------------------------------------------------------------------------------- */
void setlog(const char fname[])
{
    debuglog("Init[LOG]:\t\t Starting log to %30s \t\t\t\t\t\t\t\t\t\t\t\t\t\t\t[done]\n", fname);
	remove(fname);
    rename(_wlog_file, fname);
    strncpy(_wlog_file, fname, FILENAME_LENGTH - 1);
}
/* ------------------------------------------------------------------------------- */
void setlogsilent(const char fname[])
{
	remove(fname);
	strncpy(_wlog_file, fname, FILENAME_LENGTH - 1);
}
/* ------------------------------------------------------------------------------------------ */
void marklog()
{
	writelog("=============================================================================================================================\n");
}
/* ------------------------------------------------------------------------------------- */
void marklog2()
{
	writelog("-----------------------------------------------------------------------------------------------------------------------------\n");
}
/* ------------------------------------------------------------------------------------- */
int MeasureTime(int Reset)
{
    static struct timeval t1;
    static struct timeval t3;
    static struct timeval t5;
    static struct timeval t7;
    struct timeval        tpom;
    int tsec, tmics;
    tmics = tsec = 0;
    switch (Reset)
    {
    default:
    case 0: gettimeofday(&t1, NULL); break; /* LINEAR SOLVER */
    case 2: gettimeofday(&t3, NULL); break; /* NON-LINEAR SOLVER */
    case 4: gettimeofday(&t5, NULL); break; /* DISCRETIZATION */
    case 6: gettimeofday(&t7, NULL); break; /* OTHER */
    case 1:
        gettimeofday(&tpom, NULL);
        tsec = tpom.tv_sec - t1.tv_sec;
        if (tpom.tv_usec > t1.tv_usec) tmics = tpom.tv_usec - t1.tv_usec;
        else 
        { 
	    tsec--;  tmics = 1000000 + tpom.tv_usec - t1.tv_usec;
        }
        writelog("[%d:%02d.%02ds]",tsec / 60, tsec % 60, tmics / 10000);
        break;
    case 3:
        gettimeofday(&tpom,NULL);
        tsec = tpom.tv_sec - t3.tv_sec;
        if (tpom.tv_usec > t3.tv_usec)
            tmics = tpom.tv_usec - t3.tv_usec;
        else
        {
            tsec--;
            tmics = 1000000 + tpom.tv_usec - t3.tv_usec;
	} 
        writelog("[%d:%02d.%02ds]",tsec / 60, tsec % 60, tmics / 10000);
        break;
    case 5:
        gettimeofday(&tpom,NULL);
        tsec = tpom.tv_sec - t5.tv_sec;
        if (tpom.tv_usec > t5.tv_usec)
	    tmics = tpom.tv_usec - t5.tv_usec;
	else
	{
	    tsec--;
	    tmics = 1000000 + tpom.tv_usec - t5.tv_usec;
	} 
        writelog("[%d:%02d.%02ds]",tsec / 60, tsec % 60, tmics / 10000);
        break;
    case 7:
        gettimeofday(&tpom,NULL);
        tsec = tpom.tv_sec - t7.tv_sec;
        if (tpom.tv_usec > t7.tv_usec)
            tmics = tpom.tv_usec - t7.tv_usec;
        else
        {
            tsec--;
            tmics = 1000000 + tpom.tv_usec - t7.tv_usec;
        } 
        writelog("[%d:%02d.%02ds]",tsec / 60, tsec % 60, tmics / 10000);
        break;
    }
    return tsec * 100 + tmics / 10000;
}
/* ------------------------------------------------------------------------------- */
void MemCheck(void *a, const char format[],...)
{
    va_list argptr;
    FILE    *log;

    if (a == NULL) 
    { 

        writelog("\nError[MEMORY]: Critical error, memory not allocated!\n");  
        log = fopen(_wlog_file, "a+"); 
        va_start(argptr, format);
        vfprintf(log, format, argptr);  
        va_end(argptr);
        va_start(argptr, format);
        vfprintf(stdout, format, argptr);
        va_end(argptr);
        fclose(log);
        exit(1); 
    }
}
/* ------------------------------------------------------------------------------- */
void FileCheck(void *a, const char format[],...)
{
    va_list argptr;
    FILE    *log;
    if (a == NULL) 
    { 
        writelog("\nError[FILE]: Critical error, file opening failed!\n"); 
        log = fopen(_wlog_file, "a+"); 
        va_start(argptr, format);
        vfprintf(log, format, argptr);  
        va_end(argptr);
        va_start(argptr, format);
        vfprintf(stdout, format, argptr);
        va_end(argptr);
        fclose(log);
        exit(1);
    }
}
/* ------------------------------------------------------------------------------- */
void CriticalErrorMsg(const char format[],...)
{ 
    FILE    *log;
    va_list  argptr;
    if (strlen(_wlog_file) == 0)
        setlogsilent("writelog.log");
    log = fopen(_wlog_file, "a+"); 
    if (log != NULL)
    {
        va_start(argptr, format);
        vfprintf(log, format, argptr);  
        va_end(argptr);
        va_start(argptr, format);
        vfprintf(stdout, format, argptr);
        va_end(argptr);
        fclose(log);
    }
    exit(1); 
}
/* ----------------------------------------------------------------------------------------- */
void debuglog(const char format[],...)
{
#ifdef DEBUG
    FILE    *log;
    va_list  argptr;
    if (strlen(_wlog_file) == 0)
        setlogsilent("writelog.log");
    log = fopen(_wlog_file, "a+"); 
    if (log != NULL)
    {
        va_start(argptr, format);
        vfprintf(log, format, argptr);  
        va_end(argptr);
        va_start(argptr, format);
        vfprintf(stdout, format, argptr);
        va_end(argptr);
        fclose(log);
    }
#endif
}
/* ----------------------------------------------------------------------------------------- */
