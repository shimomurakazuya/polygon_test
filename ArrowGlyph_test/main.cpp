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
#include "ArrowGlyph.h"
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
    int n_glyphs = 2;
    // 4 debug
    kvs::ValueArray<kvs::Real32> position;
    position.allocate(3*n_glyphs);
    position.at(0) = 0;
    position.at(1) = 0;
    position.at(2) = 0;
    position.at(3) = 2;
    position.at(4) = 2;
    position.at(5) = 2;
    kvs::ValueArray<kvs::Real32> direction;
    direction.allocate(3*n_glyphs);
    direction.at(0) = 1;
    direction.at(1) = 0;
    direction.at(2) = 0;
    direction.at(3) = -1;
    direction.at(4) = 0;
    direction.at(5) = 0;

    // Create an sphere glyph polygon.
//    kvs::PolygonObject* glyph_polygon = new kvs::SphereGlyph( n_glyphs, position );

      kvs::PolygonObject* glyph_polygon = new kvs::ArrowGlyph( n_glyphs, position, direction );
    kvs::glut::Screen screen( &app );
    screen.registerObject( glyph_polygon );
    screen.show();

    return app.run();
}
