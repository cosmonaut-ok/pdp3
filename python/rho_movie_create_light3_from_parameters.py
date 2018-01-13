#!/usr/bin/env python3

import os
import argparse

from xml.dom import minidom
from numpy import *
from pylab import *

# from matplotlib import *
import matplotlib.pyplot as plt
import matplotlib.animation as ani

# import matplotlib.pyplot as pyplot
# import matplotlib.animation as animation

class Pdp3Movie:
    def __init__(self, xml_config_file, video_file=None, file_delta=100, clim_e1=[0,1], clim_e3=[0,1], clim_rho_beam=[-1e-7, 0]):
        self.xml_config_file = xml_config_file
        self.file_delta = file_delta
        self.clim_e1 = clim_e1
        self.clim_e3 = clim_e3
        self.clim_rho_beam = clim_rho_beam
        
        self.__config_path = os.path.dirname(xml_config_file)
        # calculate movie path
        video_path = os.path.dirname(video_file) if video_file else self.__config_path
        video_file_name = os.path.basename(video_file) if video_file else 'field_movie.avi'
        self.__movie_filename = os.path.join(video_path, video_file_name)

        ## read xml_config_file
        dom_root = minidom.parse(xml_config_file)

        # get geometry parameters for gird and ticks
        geometry = dom_root.getElementsByTagName('geometry')[0]
        self.__size_1 = int(geometry.getElementsByTagName('n_grid_r')[0].firstChild.data)-1
        self.__size_3 = int(geometry.getElementsByTagName('n_grid_z')[0].firstChild.data)-1
        y_tick_max = float(geometry.getElementsByTagName('r_size')[0].firstChild.data)
        x_tick_max = float(geometry.getElementsByTagName('z_size')[0].firstChild.data)

        # %% set number of marks with names in X and Y axes
        x_ticks_number = 10;
        y_ticks_number = 4;
        # %% Titles and Ticks
        self.__x_axe_title = 'Z(m)';
        self.__y_axe_title = 'R(m)';
        self.__x_tick_range = linspace(0, x_tick_max, x_ticks_number) # we need 10 (or x_ticks_number) ticks
        self.__x_tick_gird_size = linspace(0, self.__size_3, x_ticks_number) # from 0 to x_tick_max. it's required 
        self.__y_tick_range = linspace(0, y_tick_max, y_ticks_number) # to convert gird to real size (meters)
        self.__y_tick_gird_size = linspace(0, self.__size_1, y_ticks_number) # Same for X and Y axes
        
        # get file to save parameters
        file_save_parameters = dom_root.getElementsByTagName('file_save_parameters')[0]
        local_data_path = file_save_parameters.getElementsByTagName('path_to_result')[0].firstChild.data

        # calculate data path
        if str.startswith(local_data_path, '/'):
            self.__data_path = local_data_path
        else:
            self.__data_path = os.path.join(self.__config_path, local_data_path)

        self.__data_file_e1 = os.path.join(self.__data_path, 'e1')
        self.__data_file_e3 = os.path.join(self.__data_path, 'e3')
        self.__data_file_rho_beam = os.path.join(self.__data_path, 'rho_beam')
        self.__figure = plt.figure()

    def __make_figure(self):
        cmap = 'gray'
        
        fig = plt.figure()
        fig.set_size_inches([10,5])
        
        tight_layout()
        
        a1 = fig.add_subplot(311)
        a1.set_aspect('equal')
        ax1 = a1.get_xaxis()
        ay1 = a1.get_yaxis()
        im1 = a1.imshow(rand(self.__size_1,self.__size_3),cmap=cmap,interpolation='nearest')
        im.set_clim(self.clim_e1)
        
        a2 = fig.add_subplot(312)
        a2.set_aspect('equal')
        ax2 = a2.get_xaxis()
        ay2 = a2.get_yaxis()
        im2 = a2.imshow(rand(self.__size_1,self.__size_3),cmap=cmap,interpolation='nearest')
        im.set_clim(self.clim_e3)
        
        a3 = fig.add_subplot(313)
        a3.set_aspect('equal')
        ax3 = a3.get_xaxis()
        ay3 = a3.get_yaxis()
        im3 = a3.imshow(rand(self.__size_1,self.__size_3),cmap=cmap,interpolation='nearest')
        im.set_clim(self.clim_rho_beam)

        return fig, im1, im2, im3, a1, a2, a3

    def __aniframe(self, data1, data2, data3):
        fig, im1, im2, im3 = self.__make_figure()
        # ax = fig.add_subplot(111)
        # ax.set_aspect('equal')
        # ax.get_xaxis().set_visible(False)
        # ax.get_yaxis().set_visible(False)
        
        # im = ax.imshow(rand(self.__size_1,self.__size_3),cmap='gray',interpolation='nearest')
        
        # im.set_clim([0,1])
        
        def update_img(n):
            local_step = n % self.file_delta ## mod(t, file_delta)
            print("Processing frame %d" % (local_step))
            
            h_field_shot1 = data[self.__size_1*self.__size_3*local_step:self.__size_1*self.__size_3*(local_step+1)]
            h_field_matrix1 = fliplr(reshape(h_field_shot1, (self.__size_1,self.__size_3)))
            
            h_field_shot2 = h_field_e3[self.__size_1*self.__size_3*local_step:self.__size_1*self.__size_3*(local_step+1)]
            h_field_matrix2 = fliplr(reshape(h_field_shot2, (self.__size_1,self.__size_3)))

            h_field_shot3 = data[self.__size_1*self.__size_3*local_step:self.__size_1*self.__size_3*(local_step+1)]
            h_field_matrix3 = fliplr(reshape(h_field_shot3, (self.__size_1,self.__size_3)))


            # tmp = rand(300,300)
            im1.set_data(h_field_matrix1)
            im1.set_data(h_field_matrix2)
            return im
        
        #legend(loc=0)
        anim = ani.FuncAnimation(fig,update_img,100,interval=30)
        writer = ani.writers['ffmpeg'](fps=10)
        dpi = 300
        anim.save(self.__movie_filename, writer=writer, dpi=dpi)
        return ani

    def create_movie_with_3_plots(self):
        # % create figure window (and place it as current figure)

        
        N = self.file_delta #  100 by default
        
        for k in range(0, 1): # TODO: why to 100?
            print("Loading files set %d" % (k))
            
            tstart = k*N # k=2 -> 200
            tend = ((k+1)*N-1) # k=2 -> 3*100-1 -> 299
            i = 1
            
            if not os.path.isfile(self.__data_file_e1 + str(k)) \
               or not os.path.isfile(self.__data_file_e3 + str(k)) \
               or not os.path.isfile(self.__data_file_rho_beam + str(k)):
                print('No more data files exists. Exiting')
                return
            
            ## Open data files
            fidh_e1 = open(self.__data_file_e1 + str(k), 'r')
            fidh_e3 = open(self.__data_file_e3 + str(k), 'r')
            fidh_rho_beam = open(self.__data_file_rho_beam + str(k), 'r')
            
            # %e is an exponential notation: ex. 5e2
            # h_field_e1 = fscanf(fidh_e1, '%e', size_1*size_3*file_delta); # TODO: why do we times all values?
            h_field_e1 = fromfile(fidh_e1, dtype=float, count=self.__size_1*self.__size_3*self.file_delta, sep=' ')
            # h_field_e3 = fscanf(fidh_e3, '%e', size_1*size_3*file_delta);
            ## h_field_e3 = fromfile(fidh_e3, dtype=float, count=self.__size_1*self.__size_3*self.file_delta, sep=' ')
            # h_field_rho_beam = fscanf(fidh_rho_beam, '%e', size_1*size_3*file_delta);
            ## h_field_rho_beam = fromfile(fidh_rho_beam, dtype=float, count=self.__size_1*self.__size_3*self.file_delta, sep=' ')
            ## Close data files
            ## fidh_e1.close()
            ## fidh_e3.close()
            fidh_rho_beam.close()

            self.__aniframe(h_field_e1)
            return 0

            with writer.saving(fig, self.__movie_filename, 100):
                for t in range(tstart, tend):
                    local_step = t % self.file_delta ## mod(t, file_delta)
                    
                    print("Processing frame %d" % (local_step))
                    
                    h_field_shot1 = h_field_e1[self.__size_1*self.__size_3*local_step+1:self.__size_1*self.__size_3*(local_step+1)+1]
                    h_field_matrix1 = fliplr(reshape(h_field_shot1, (self.__size_1,self.__size_3)))
                    
                    ## h_field_shot2 = h_field_e3[self.__size_1*self.__size_3*local_step+1:self.__size_1*self.__size_3*(local_step+1)+1]
                    ## h_field_matrix2 = fliplr(reshape(h_field_shot2, (self.__size_1,self.__size_3)))
                    
                    # h_field_shot3 = h_field_rho_beam[self.__size_1*self.__size_3*local_step+1:self.__size_1*self.__size_3*(local_step+1)+1]
                    # h_field_matrix3 = fliplr(reshape(h_field_shot3, (self.__size_1,self.__size_3)))
                    
                    ax.imshow(h_field_matrix1)
                    writer.grab_frame()
                    
                    #     set(im1, 'CData', h_field_matrix1);
                    #     set(im2, 'CData', h_field_matrix2);
                    #     set(im3, 'CData', h_field_matrix3);
                    
                    #     figHandle = gcf;
                    
                    #     frame = getframe(figHandle);
                    #     D(:,:,:,i)   = frame.cdata;
                    #     D(:,:,:,i+1) = frame.cdata;
                    #     D(:,:,:,i+2) = frame.cdata;
                    #     i = i + 3;
                    #     drawnow;
                    
                    # f2 = figure('Position', [20 20 100 100]);
                    # at = axes('Parent', f2);
                    # mov = immovie(D);
                    # %% movie_filename = strcat(video_path,'field_movie_',num2str(t/1e-7,'%3.2f'),'.avi');
                    # %% movie_object = VideoWriter(movie_filename);
                    # %% open(movie_object);
                    # writeVideo(movie_object, mov);
                    # clear mov;
                    # %% close(movie_object)
                    # close(f2);
                    # %% pack


def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('properties_path', metavar='properties_path', type=str,
                        help='Full path to properties.xml')

    parser.add_argument('--video_file', type=str,
                        help='Full path to output video file')

    # parser.add_argument('clim_e1', metavar='properties_path', type=str,
    #                     help='Full path to properties.xml')

    ## video_file=None, file_delta=100, clim_e1=[0,1], clim_e3=[0,1], clim_rho_beam=[-1e-7, 0]):


    parser.add_argument('--clim_e1', type=str,
                        help='Color limit range for e1. Default: 0:1. Not implemented')
    parser.add_argument('--clim_e3', type=str,
                        help='Color limit range for e3. Default: 0:1. Not implemented')
    parser.add_argument('--clim_rho_beam', type=str,
                        help='Color limit range for rho_beam. Default: -1e-7:0. Not implemented')

    args = parser.parse_args()
    print(args)
    
    clim_e1 = list(map(int, args.clim_e1.split(':'))) if args.clim_e1 else [0,1]
    clim_e3 = list(map(int, args.clim_e3.split(':'))) if args.clim_e3 else [0,1]
    clim_rho_beam = list(map(int, args.clim_rho_beam.split(':'))) if args.clim_rho_beam else [-1e-7, 0]
    file_delta = 100

    movie = Pdp3Movie(args.properties_path, video_file=args.video_file, file_delta=file_delta, clim_e1=clim_e1, clim_e3=clim_e3, clim_rho_beam=clim_rho_beam)
    movie.create_movie_with_3_plots()


        
if __name__ == "__main__":
    main()
