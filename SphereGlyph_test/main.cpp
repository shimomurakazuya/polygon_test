/*****************************************************************************/
/**
 *  @file   main.cpp
 *  @brief  Example program for kvs::SphereGlyph class.
 *  @author Naohisa Sakamoto
 */
/*****************************************************************************/
#include <kvs/Message>
#include <kvs/StructuredVolumeObject>
#include <kvs/StructuredVolumeImporter>
//#include <kvs/SphereGlyph>
#include <kvs/TornadoVolumeData>
#include <kvs/glut/Application>
#include <kvs/glut/Screen>
#include "SphereGlyph.h"
//#include "SphereGlyph_test.h"

/*===========================================================================*/
/**
 *  @brief  Main function.
 *  @param  argc [i] argument counter
 *  @param  argv [i] argument values
 *  @return true, if the main process is done succesfully
 */
/*===========================================================================*/
int main( int argc, char** argv )
{
    kvs::glut::Application app( argc, argv );

    /* Read volume data from the specified data file. If the data file is not
     * specified, tornado volume data is created by using kvs::TornadoVolumeData class.
     */

    // グリフポリゴン表示に関係あるのは glyphのコンストラクターとsetpolygon()と各種セッターのみ
    int n_glyphs = 1;
    // 4 debug
    kvs::ValueArray<kvs::Real32> pos;
    pos.allocate(3*n_glyphs);
    pos.at(0) = 0;
    pos.at(1) = 0;
    pos.at(2) = 0;
    //pos.at(3) = 0;
    //pos.at(4) = 0;
    //pos.at(5) = 3;
  
    // Create an sphere glyph polygon.
    kvs::PolygonObject* glyph_polygon = new kvs::SphereGlyph( n_glyphs, pos );
    kvs::glut::Screen screen( &app );
    screen.registerObject( glyph_polygon );
    screen.show();

    return app.run();
}
