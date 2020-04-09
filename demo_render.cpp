/*
  ======================================================================
   demo.cpp - SFML-enabled version of the demo that provides visual feedback
  ----------------------------------------------------------------------
   Author : Sean Mcfarlane
   Creation Date : Jan 2019

   Description:

    This project was created based on the paper "Real-Time Fluid Dynamics for Games" by Jos Stam.
    Code is a modified version of the sample code provided alongside the paper, which was adapted to use C++ and SFML rather than C and OpenGL.
    This was created for the purpose of testing optimization methods as part of my CSC485C coursework.

  =======================================================================
*/

#include <stdlib.h>
#include <stdio.h>
#include <SFML/Graphics.hpp>
#include "project.h"
#include "cudatest.h"

/* global variables */

static uint32_t dvel;

uint8_t *pixels;

static int win_x, win_y;
static uint32_t mouse_down[3];
static int omx, omy, mx, my;

sf::RenderWindow *window;
sf::Clock *g_clock;
sf::Font *font;

bool CKeyDown;
uint32_t render_mode;
int brush_size;

const float visLen = 30.0f;

/*
  ----------------------------------------------------------------------
   SFML specific routines
  ----------------------------------------------------------------------
*/

static void init_sfml(uint32_t resolution)
{	
	font = new sf::Font();
	font->loadFromFile("UbuntuMono-B.ttf");
	sf::ContextSettings settings;
	settings.antialiasingLevel = 4;
	window = new sf::RenderWindow(sf::VideoMode(resolution, resolution), "Real-Time Fluid Dynamics Sim", sf::Style::Default, settings);
	window->setFramerateLimit(frameRate);
	g_clock = new sf::Clock();
	pixels = new uint8_t[resolution*resolution*4];
}

static void mouse_func(uint32_t id, bool pressed)
{
	mouse_down[id] = pressed;
	
	sf::Vector2i mPos = sf::Mouse::getPosition(*window);

	omx = mx = mPos.x;
	omx = my = mPos.y;

}

static void motion_func()
{
	sf::Vector2i mPos = sf::Mouse::getPosition(*window);
	mx = mPos.x;
	my = mPos.y;
}

static void render_density ( void )
{
	uint32_t texsize = resolution;

	sf::Texture tex;
	if (!tex.create(texsize, texsize)) std::cout << "tex broked";
	sf::Sprite sprite(tex);

	for (int x = 0; x < texsize; x++)
	{
		for (int y = 0; y < texsize; y++)
		{
			float xStep = (float)x / (float)texsize;
			float yStep = (float)y / (float)texsize;

			uint32_t gridx = xStep * N + pad;
			uint32_t gridy = (N-(yStep * N));

			float density = dens_prev[IX(gridx, gridy)];

			density = (density > 255) ? 255 : density; //Clamp to 255.

			pixels[((y * texsize) + x) * 4 + 0] = (uint8_t)255-density;	// R
			pixels[((y * texsize) + x) * 4 + 1] = (uint8_t)255-density;	// G
			pixels[((y * texsize) + x) * 4 + 2] = (uint8_t)255-density;	// B
			pixels[((y * texsize) + x) * 4 + 3] = (uint8_t)255;			// A
		}
	}
	tex.update(pixels);

	window->draw(sprite);
	window->display();
}


static void Render(uint32_t optim_mode){
	if (render_mode == 0) {
		switch (optim_mode) {
			case 0: // Unoptimized
			{
				base::render_velocity();
				break;
			}
			case 1: // Single core optimized (ZIX)
			{
				opt::render_velocity();
				break;
			}
			case 2: // Parallelized (ZIX)
			{
				opt::render_velocity();
				break;
			}
			case 3: // SIMD (IX)
			{
				SIMD::render_velocity();
				break;
			}
			case 4: // SIMD & Parallelized (IX)
			{
				SIMD::render_velocity();
				break;
			}
			case 5: // CUDA
			{
				CUDA::render_velocity();
				break;
			}
		}
	}
	else if (render_mode == 1) {
		render_density();
	}
	else {
		render_mode = 0;
		//std::cout << "Set render mode back to 0\n";
	}
}

static void ReadInput(float *d, float *u, float *v)
{
	uint32_t i, j, size = array_size;

	if (sf::Keyboard::isKeyPressed(sf::Keyboard::C)) 
	{
		if (!CKeyDown) 
		{
			render_mode++;
			CKeyDown = true;
		}
	}
	else if(CKeyDown)
	{
		CKeyDown = false;
	}

	for (i = 0; i < size; i++)
	{
		u[i] = v[i] = d[i] = 0.0f;
	}

	if (!mouse_down[0] && !mouse_down[2])
		return;

	i = (int)((mx / (float)win_x) * N)+pad;
	j = (int)(((win_y - my) / (float)win_y) * N)+pad;

	if (i < pad || i > N || j < pad || j > N)
		return;


	if (mouse_down[0]){
		if(optim_mode==1 || optim_mode==2){
			if (brush_size > 1) {
				for (int i_offset = -brush_size / 2; i_offset <= brush_size / 2; i_offset++)
				{
					for (int j_offset = -brush_size / 2; j_offset <= brush_size / 2; j_offset++)
					{
						if (i_offset * i_offset + j_offset * j_offset > (brush_size / 2) * (brush_size / 2)) { continue; }
						uint32_t index = ZIX(i + i_offset, j + j_offset);
						if (index < array_size && index >= 0) {
							u[index] = force * dt * (mx - omx);
							v[index] = force * dt * (omy - my);
						}
					}
				}
			}
		}
		else{

			if (brush_size > 1) {
				for (int i_offset = -brush_size / 2; i_offset <= brush_size / 2; i_offset++){
					for (int j_offset = -brush_size / 2; j_offset <= brush_size / 2; j_offset++){
						if (i_offset * i_offset + j_offset * j_offset > (brush_size / 2) * (brush_size / 2)) { continue; }
						uint32_t index = IX(i + pad + i_offset, j + pad + j_offset);
						if (index < array_size && index >= 0) {
							u[index] = force * dt * (mx - omx);
							v[index] = force * dt * (omy - my);
						}
					}
				}
			}
		}
	}

	if (mouse_down[2]){
		if(optim_mode==1 || optim_mode==2){
			if (brush_size > 1) {
				for (int i_offset = -brush_size / 2; i_offset <= brush_size / 2; i_offset++){
					for (int j_offset = -brush_size / 2; j_offset <= brush_size / 2; j_offset++){
						if (i_offset * i_offset + j_offset * j_offset > (brush_size / 2) * (brush_size / 2)) { continue; }
						uint32_t index = ZIX(i + i_offset, j + j_offset);
						if (index < array_size && index >= 0) {
							d[index] += source;
						}
					}
				}
			}
		}
		else{
			if (brush_size > 1) {
				for (int i_offset = -brush_size / 2; i_offset <= brush_size / 2; i_offset++){
					for (int j_offset = -brush_size / 2; j_offset <= brush_size / 2; j_offset++){
						if (i_offset * i_offset + j_offset * j_offset > (brush_size/2) * (brush_size/2)) { continue; }
						uint32_t index = IX(i + pad + i_offset, j + pad + j_offset);
						if (index < array_size && index >= 0) {
							d[index] += source;
						}
					}
				}
			}
		}
	}

	omx = mx;
	omy = my;

	return;
}

static void GameLoop(){
	DPRINT("GameLoop begin\n");
	while (window->isOpen())
	{
		cur_iter++;
        DPRINT("Iteration "<< cur_iter <<"\n");
		sf::Event event;
		while (window->pollEvent(event))
		{
			if (event.type == sf::Event::Closed){
				DPRINT("Closed\n");
				window->close();
			}
			if (event.type == sf::Event::MouseButtonPressed){
				if(event.mouseButton.button == sf::Mouse::Left){
					mouse_func(0, true);
				}
				if (event.mouseButton.button == sf::Mouse::Middle){
					mouse_func(1, true);
				}
				if (event.mouseButton.button == sf::Mouse::Right){
					mouse_func(2, true);
				}
			}
			if (event.type == sf::Event::MouseButtonReleased){
				if (event.mouseButton.button == sf::Mouse::Left){
					mouse_func(0, false);
				}
				if (event.mouseButton.button == sf::Mouse::Middle){
					mouse_func(1, false);
				}
				if (event.mouseButton.button == sf::Mouse::Right){
					mouse_func(2, false);
				}
			}
			if (event.type == sf::Event::MouseWheelMoved)
			{
				brush_size += event.mouseWheel.delta;
			}
			if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::R)
			{
				clear_data();
			}
		}

		motion_func();

		ReadInput(dens_prev, u_prev, v_prev);
		Simulate(optim_mode);
		Render(optim_mode);
	}
}

namespace opt
{
	void render_velocity(){
		DPRINT("render_velocity begin\n");
		sf::VertexArray line(sf::Lines, 2 * (N + 2) * (N + 2));

		DPRINT("1\n");
		#ifdef DISPLAY_COORDS
		sf::Text texts[(N+bnd)*(N+bnd)];
		#endif
		DPRINT("2\n");
		uint32_t hits = 0;
		uint32_t zX, zY, i, j;

		DPRINT("RENDERLOOP: zonesInRow = (" << zonesInRow <<")\n");
		DPRINT("RENDERLOOP: zoneLen = (" << zoneLen <<")\n");
		DPRINT("RENDERLOOP: ARRAY SIZE = (" << 2*N*N <<")\n");


		for( zY=0; zY<zonesInRow; zY++){
			for( zX=0; zX<zonesInRow; zX++){
				for ( j=zY*zoneLen; j<(zY+1)*zoneLen; j++ ) {
					for ( i=zX*zoneLen; i<(zX+1)*zoneLen; i++ ) { 
						if(i <= 0 || j <= 0 || i >= N || j >= N){
							DPRINT("RENDERLOOP: (" << i << "," << j <<") OMITTED\n");	
							#ifdef DISPLAY_COORDS				
							uint32_t zix = ZIX(i, j);
							float xStep = (float)resolution/N;
							float yStep = (float)resolution/N;	
							texts[zix] = sf::Text(sf::String("("+std::to_string(i)+","+std::to_string(j)+")"), *font);
							texts[zix].setCharacterSize(11);
							texts[zix].setFillColor(sf::Color::Red);
							sf::FloatRect textRect = texts[zix].getLocalBounds();
							float xPos = xStep * ((float)i);
							float yPos = resolution-(yStep * ((float)j));
							xPos -= textRect.width/2;
							yPos -= textRect.height/2;
							texts[zix].setPosition(xPos, yPos);
							#endif
							continue;
						}

                
						uint32_t zix = ZIX(i, j);
						uint32_t index = 2 * zix;
						DPRINT("RENDERLOOP: (" << i << "," << j <<")="<< index <<"\n");

						float xStep = (float)resolution/N;
						float yStep = (float)resolution/N;

						float xPos = xStep * ((float)i);
						float yPos = resolution-(yStep * ((float)j));
						float xOffset = u[zix]*visLen;
						float yOffset = v[zix]*visLen;

						float colInterp = (abs(xOffset)+abs(yOffset))/(visLen*0.5f);
						colInterp = colInterp <= 1 ? colInterp : 1.0f;

						line[index].position = sf::Vector2f(xPos, yPos);
						line[index].color.r = 255*colInterp;
						line[index].color.b = 255*(1-colInterp);
						line[index].color.g = 128;
						line[index+1].position = sf::Vector2f(xPos + xOffset,  yPos + yOffset);
						line[index+1].color = line[index].color;
						#ifdef DISPLAY_COORDS
						texts[zix] = sf::Text(sf::String("("+std::to_string(i)+","+std::to_string(j)+")"), *font);
						texts[zix].setCharacterSize(11);
						texts[zix].setFillColor(line[index].color);
						sf::FloatRect textRect = texts[zix].getLocalBounds();
						texts[zix].setPosition(xPos - textRect.width/2, yPos- textRect.height/2);
						#endif
					}
				}	
			}
		}

		window->clear();
		DPRINT("Cleared window \n");
		#ifdef DISPLAY_COORDS
		if(DISPLAY_COORDS){
			for(uint32_t i = 0; i < (N+2)*(N+2); i++){
				window->draw(texts[i]);
			}
		}
		#endif
		window->draw(line);
		window->display();
	}
} // namespace opt

namespace base
{
	void render_velocity(){
		sf::VertexArray line(sf::Lines, 2 * N * N);
		for(uint32_t y = 0; y < N; y++){
			for(uint32_t x = 0; x < N; x++){
				uint32_t index = 2*(x+(y*N));
				//DPRINT("("<< x*2 << ", " << y << ")=ArrayPos " << index << "\n"); 
				//DPRINT("("<< x*2+1 << ", " << y << ")=ArrayPos " << index+1 << "\n");
				
				float xStep = (float)resolution/N;
				float yStep = (float)resolution/N;

				float xPos = xStep * ((float)x + 0.5f);
				float yPos = yStep * ((float)y + 0.5f);		

				float xOffset = u[IX(x, N-y)]*visLen;
				float yOffset = v[IX(x, N-y)]*visLen;
				float colInterp = (abs(xOffset)+abs(yOffset))/(visLen*0.5f);
				colInterp = colInterp <= 1 ? colInterp : 1.0f;

				line[index].color.r = 255*colInterp;
				line[index].color.b = 255*(1-colInterp);
				line[index].color.g = 128;
				
				line[index].position = sf::Vector2f(xPos, yPos);

				line[index+1].position = sf::Vector2f(xPos+xOffset, yPos+yOffset);
				line[index+1].color = line[index].color;
			}
		}
		window->clear();
		window->draw(line);
		window->display();
	}
} // namespace base


namespace SIMD
{
	void render_velocity(){
		sf::VertexArray line(sf::Lines, 2 * N * N);
		for(uint32_t y = 0; y < N; y++){
			for(uint32_t x = 0; x < N; x++){
				uint32_t index = 2*(x+(y*N));
				//DPRINT("("<< x*2 << ", " << y << ")=ArrayPos " << index << "\n"); 
				//DPRINT("("<< x*2+1 << ", " << y << ")=ArrayPos " << index+1 << "\n");
				
				float xStep = (float)resolution/N;
				float yStep = (float)resolution/N;

				float xPos = xStep * ((float)x + 0.5f);
				float yPos = yStep * ((float)y + 0.5f);		

				float xOffset = u[IX(x+pad, N+pad-y)]*visLen;
				float yOffset = v[IX(x+pad, N+pad-y)]*visLen;
				float colInterp = (abs(xOffset)+abs(yOffset))/(visLen*0.5f);
				colInterp = colInterp <= 1 ? colInterp : 1.0f;

				line[index].color.r = 255*colInterp;
				line[index].color.b = 255*(1-colInterp);
				line[index].color.g = 128;
				
				line[index].position = sf::Vector2f(xPos, yPos);

				line[index+1].position = sf::Vector2f(xPos+xOffset, yPos+yOffset);
				line[index+1].color = line[index].color;
			}
		}
		window->clear();
		window->draw(line);
		window->display();
	}
} // namespace SIMD

namespace CUDA
{
	void render_velocity() {
		base::render_velocity();
	}
} // namespace CUDA


/*
  ----------------------------------------------------------------------
   main --- main routine
  ----------------------------------------------------------------------
*/

int main ( int argc, char ** argv )
{

	if ( argc != 1 && argc != 4 ) {
		fprintf ( stderr, "Usage: demo.exe <select implementation[0-5]> <resolution[int]> <iterations[int]>\n");
		return 1;
	}

	profiling = false;
	if ( argc == 1 ) {
		optim_mode = 0;
		N = 256;
		iterations = 1000;
		fprintf ( stderr, "Using defaults: profiler mode, SIMD, 256x256, 1000 iterations\n");
	} else {
		optim_mode = atoi(argv[1]);
		N = atoi(argv[2]);
		iterations = atoi(argv[3]);
	}

	if(optim_mode==3||optim_mode==4){bnd = 8;}
	else{bnd = 2;}
	N-=bnd;
    pad = bnd/2;
	array_size = (N + bnd) * (N + bnd);
	timeSpeed = 1.0f;
	diff = 0.00001f;
	visc = 0.00001f;
	force = 50.0f;
	source = 512.0f;
	brush_size = resolution/128;

	printf ( "\n\nHow to use this demo:\n\n" );
	printf ( "\t Add densities with the right mouse button\n" );
	printf ( "\t Add velocities with the left mouse button and dragging the mouse\n" );
	printf ( "\t Toggle density/velocity display with the 'c' key\n");
	printf ( "\t Use the scroll wheel to adjust brush size.\n" );
	printf ( "\t Clear the simulation by pressing the 'r' key\n" );
	printf ( "\t Quit by closing the window or using CTRL+C on the command line.\n" );

	dvel = 0;

	if (optim_mode == 5) {
		if (!allocate_data()) return 1;
	}
	else {
		if (!allocate_data_simd()) return 1;
	}

	clear_data();

	win_x = resolution;
	win_y = resolution;

	dt = timeSpeed / 60.0f;

	init_sfml(resolution);
	GameLoop();

	return 0;
}