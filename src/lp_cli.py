#!python3

try:
	import argparse
	import local_pcangsd as lp
except ImportError:
	raise ImportError('Something went wrong in importing modules')


def main():

	parser = argparse.ArgumentParser(prog="local_pcangsd", description="Runs local_pcangsd")
	# parser.add_argument()

	print('hello world')


if __name__ == "__main__":
	main()