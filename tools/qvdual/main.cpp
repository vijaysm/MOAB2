#include <qapplication.h>
#include "uiQVDual.h"
#include "MBInterface.hpp"

MBInterface *gMB = NULL;

int main( int argc, char ** argv )
{
    QApplication a( argc, argv );
    uiQVDual w;
    w.show();
    a.connect( &a, SIGNAL( lastWindowClosed() ), &a, SLOT( quit() ) );

    if (argc > 1) w.fileOpen(QString(argv[1]));
    
    return a.exec();
}
