#ifndef BASILISK_HEADER_11
#define BASILISK_HEADER_11
#line 1 "/home/dahuanghhc/basilisk/src/gl/font.h"
/*
 * Portions copyright (C) 2004, the OpenGLUT project contributors.
 * OpenGLUT branched from freeglut in February, 2004.
 *
 * Copyright (c) 1999-2000 Pawel W. Olszta. All Rights Reserved.
 * Written by Pawel W. Olszta, <olszta@sourceforge.net>
 * Creation date: Thu Dec 16 1999
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * PAWEL W. OLSZTA BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include <stdio.h>
#include "og_font.h"

#define  GL_STROKE_ROMAN               ((void *)0x0000) // not used
#define  GL_STROKE_MONO_ROMAN          ((void *)0x0001)

/* -- IMPORT DECLARATIONS -------------------------------------------------- */

/*
 * These are the font faces defined in og_font_data.c file:
 */
/* extern SOG_StrokeFont ogStrokeRoman; */
extern SOG_StrokeFont ogStrokeMonoRoman;


/* -- PRIVATE FUNCTIONS ---------------------------------------------------- */

/*!
    Matches a font ID with a SOG_StrokeFont structure pointer.
    This was changed to match the GLUT header style.
*/
static SOG_StrokeFont *oghStrokeByID( void *font )
{
  /*    if( font == GL_STROKE_ROMAN      )
        return &ogStrokeRoman; */
    if( font == GL_STROKE_MONO_ROMAN )
        return &ogStrokeMonoRoman;

    fprintf (stderr,  "stroke font %p not found", font );
    return 0;
}


/* -- INTERFACE FUNCTIONS -------------------------------------------------- */

/*!
    \fn
    \brief    Draw a stroked character.
    \ingroup  fonts
    \param    fontID    A GLUT stroked font identifier.
    \param    character An ASCII character other than NUL.

              This function draws one \a character from one stroked font
              (selected by \a fontID)
              using OpenGL \a GL_LINE_STRIP.  These characters
              are drawn at the origin in model space.
              The the model space origin is translated at the end,
              according to the \a character width in \a fontID.

              Does nothing if:
               - The \a fontID is invalid.
               - The \a character is out of the font's range.

    \see      glBegin(), glTranslatef(), glutStrokeWidth(), glutStrokeString(),
              glutStrokeHeight(), glutBitmapCharacter()
*/
void gl_StrokeCharacter( int character )
{
    void *fontID = GL_STROKE_MONO_ROMAN;
    const SOG_StrokeChar *schar;
    const SOG_StrokeStrip *strip;
    int i, j;
    SOG_StrokeFont *font = oghStrokeByID( fontID );

    if( !font ||
        ( 1 > character ) ||
        ( font->Quantity < character ) )
        return;

    schar = font->Characters[ character ];
    if( schar )
    {
        strip = schar->Strips;

        for( i = 0; i < schar->Number; i++, strip++ )
        {
            glBegin( GL_LINE_STRIP );
            for( j = 0; j < strip->Number; j++ )
                glVertex2f( strip->Vertices[ j ].X, strip->Vertices[ j ].Y );
            glEnd( );
        }
        glTranslatef( schar->Right, 0.0, 0.0 );
    }
}

/*!
    \fn
    \brief    Draw a string of stroked characters.
    \ingroup  fonts
    \param    fontID    A GLUT stroked font identifier.
    \param    string    A NUL-terminated ASCII string.

              This function draws a \a string in the font indicated
              by \a fontID.
              It is <i>almost</i> equivalent to calling glutStrokeCharacter()
              on each character in the \a string, successively.
              Mostly, it is a convenience function to hide the loop,
              and to treat \\n as a special symbol rather than a normal
              glyph.

              The first character displays at the current model space
              origin,  The origin changes by successive translations.

              The newline character, \\n (ASCII LF) is treated as
              a newline and resets the origin horizontally
              while advancing the line 1 font-height down the y-axis.

              Does nothing if:
               - \a fontID is out of range.
               - \a string is \a NULL
               - \a string is empty

              Unlike glutBitmapString(), there is little performance
              advantage to using glutStrokeString() as compared with
              calling glutStrokeCharacter() yourself for every
              character.

    \see      glutStrokeLength(), glutStrokeCharacter(),
              glutStrokeHeight(), glutBitmapString()
*/
void gl_StrokeString( const char *string )
{
    void *fontID = GL_STROKE_MONO_ROMAN;
    int i, j;
    float length = 0.0;
    SOG_StrokeFont *font = oghStrokeByID( fontID );
    unsigned char c;

    if( font && string )
        /*
         * Step through the string, drawing each character.
         * A newline will simply translate the next character's insertion
         * point back to the start of the line and down one line.
         */
        while(( c = *string++ ))
  	    if( c < font->Quantity ) {
                if( c == '\n' )
                {
                    glTranslatef ( -length, -( float )( font->Height ), 0.0 );
                    length = 0.0;
                }
                else  /* Not an EOL, draw the bitmap character */
                {
                    const SOG_StrokeChar *schar =
                        font->Characters[ c ];
                    if( schar )
                    {
                        const SOG_StrokeStrip *strip = schar->Strips;

                        for( i = 0; i < schar->Number; i++, strip++ )
                        {
                            glBegin( GL_LINE_STRIP );
			    
                            for( j = 0; j < strip->Number; j++ )
                                glVertex2f( strip->Vertices[ j ].X,
                                            strip->Vertices[ j ].Y);

                            glEnd( );
                        }

                        length += schar->Right;
                        glTranslatef( schar->Right, 0.0, 0.0 );
                    }
                }
	    }
}

/*!
    \fn
    \brief    Returns the width in pixels of a character in a given font.
    \ingroup  fonts
    \param    fontID    A GLUT stroked font identifier.
    \param    character A character code.

              This function reports how far the model space origin will advance
              if you putput this \a character in the font named by \a fontID.
              Not all letters will use their full width, especially in
              fixed-width fonts.

              Returns 0 if \a character is out of range or if the
              \a fontID is invalid.

    \todo     Determine if any glyphs are either wider than this
              function or if they render outside of the bounding
              box given by
              <i>(0,-descent)</i> by <i>(width,height-descent)</i>.
    \note     Historically, this function has been described as
              returning a pixel-width, but was implemented to
              return the width in model-space units, rounded to integers.
              GLUT never resolved this, and freeglut duplicated the
              confusion.
              OpenGLUT has decided to stay in model-space and to
              return the unrounded floating point value.
              An unreleased GLUT 3.8 was supposed to include
              glutStrokeWidthf() and glutStrokeLengthf() (note
              the *f suffixes), but that is not in wide use.
    \see      glutStrokeCharacter(), glutStrokeLength(), glutStrokeHeight()
              glutBitmapWidth()
*/
float gl_StrokeWidth( int character )
{
    void *fontID = GL_STROKE_MONO_ROMAN;
    float ret = 0;
    SOG_StrokeFont *font = oghStrokeByID( fontID );

    if( font &&
        ( 0 < character ) &&
        ( font->Quantity > character ) )
    {
        const SOG_StrokeChar *schar = font->Characters[ character ];
        if( schar )
            ret = schar->Right;
    }

    return ret;
}

/*!
    \fn
    \brief    Returns model space width of a string in a given font.
    \ingroup  fonts
    \param    fontID    A GLUT stroked font identifier.
    \param    string    A C-style (NUL-terminated) string.

              This function reports the sum of the widths of the
              characters in a \a string, using the font metrics of
              a given \a font.

              Like glutStrokeString(), glutStrokeLength() respects
              newlines in the input.

              Returns 0 if:
               - The \a fontID is out of range.
               - The \a string is \a NULL.
               - All characters in the \a string are zero-width.

    \note     Refer to glutStrokeWidth() for notes on the
              nature of this function's return value, and for
              comparisons to old GLUT and freeglut.
    \see      glutStrokeString(), glutStrokeWidth(), glutStrokeHeight(),
              glutBitmapLength()
*/
float gl_StrokeLength( const char *string )
{
    void *fontID = GL_STROKE_MONO_ROMAN;
    unsigned char c;
    float length = 0.0;
    float this_line_length = 0.0;
    SOG_StrokeFont *font = oghStrokeByID( fontID );

    if( font && string )
        while(( c = *string++ ))
            if( c < font->Quantity )
            {
                if( c == '\n' )
                {
                    if( length < this_line_length )
                        length = this_line_length;
                    this_line_length = 0.0;
                }
                else  /* Not an EOL, increment the length of this line */
                {
                    const SOG_StrokeChar *schar =
                        font->Characters[ c ];
                    if( schar )
                        this_line_length += schar->Right;
                }
            }

    if( length < this_line_length )
        length = this_line_length;
    return length;
}

/*!
    \fn
    \brief    Returns the height of a given font.
    \ingroup  fonts
    \param    fontID    A GLUT stroked font identifier.

              This function reports the height of a font,
              given by \a fontID,
              as a global characteristic of that font.

              Returns 0 if \a fontID is invalid.

    \note     Does <i>not</i> report the height used by individual
              characters.  This may limit its usefulness; perhaps we
              should change it?  (And/or add a new function.)
    \todo     We have discussed adding a "font descender" query.
              We should go ahead and do it.
    \see      glutStrokeCharacter(), glutStrokeString(), glutStrokeWidth(),
              glutStrokeLength(), glutBitmapHeight()
*/
GLfloat gl_StrokeHeight()
{
    void *fontID = GL_STROKE_MONO_ROMAN;
    GLfloat ret = 0;
    SOG_StrokeFont *font = oghStrokeByID( fontID );

    if( font )
        ret = font->Height;

    return ret;
}

#endif
