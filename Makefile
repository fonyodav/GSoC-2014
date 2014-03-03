NAME=demo

all: $(NAME)

clean: rm $(NAME) *.o *~

run: ./$(NAME)

$(NAME): g++ $(NAME).cpp -o $(NAME) -lgivaro -llinbox -lcblas -lgmp
