mod lcgrand;
use clap::Parser;
use lcgrand::LcGrand;
use std::{cmp::Ordering, collections::BinaryHeap};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[arg(short = 'n', long = "clients", default_value = "1000")]
    clients: usize,
    #[arg(short = 'm', long = "servers", default_value = "1")]
    servers: usize,
    #[arg(short = 'b', long = "arrival")]
    arrival_mean: f64,
    #[arg(short = 'g', long = "service")]
    service_mean: f64,
    #[arg(short = 'v', long = "verbose")]
    verbose: bool,
}

#[derive(Debug)]
struct Client {
    arrival_time: f64,
    wait_time: f64,
    service_start_time: f64,
    service_time: f64,
    departure_time: f64,
}

impl Client {
    fn new(arrival_time: f64) -> Client {
        Client {
            arrival_time,
            wait_time: 0.0,
            service_start_time: 0.0,
            service_time: 0.0,
            departure_time: 0.0,
        }
    }
}

#[derive(Debug)]
enum ServerState {
    Idle,
    Busy,
}

#[derive(Debug)]
struct Server {
    state: ServerState,
    active_time: f64,
}

impl Server {
    fn new() -> Server {
        Server {
            state: ServerState::Idle,
            active_time: 0.0,
        }
    }
}

#[derive(Debug)]
enum EventState {
    Arrival,
    Departure(usize),
}

#[derive(Debug)]
struct Event {
    time: f64,
    state: EventState,
}

impl Ord for Event {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.time < other.time {
            Ordering::Greater
        } else if self.time > other.time {
            Ordering::Less
        } else {
            Ordering::Equal
        }
    }
}
impl PartialOrd for Event {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl PartialEq for Event {
    fn eq(&self, other: &Self) -> bool {
        self.time == other.time
    }
}
impl Eq for Event {}

impl Event {
    fn new(time: f64, state: EventState) -> Event {
        Event { time, state }
    }
}

struct Report {
    total_clients: usize,        // n
    total_servers: usize,        // m
    mean_interarrival_time: f64, // 1/λ
    mean_service_time: f64,      // 1/μ

    mean_wait_time: f64,
    mean_queue_length: f64,
    total_simulation_time: f64,
    servers_utilization: Vec<f64>,
}

struct Simulation {
    clock: f64,
    event_list: BinaryHeap<Event>,
    servers: Vec<Server>,
    lcgrand: LcGrand,
    arrival_mean: f64,
    service_mean: f64,
    queue: Vec<Client>,
    client_count: usize,
    served_clients: Vec<Client>,
    queue_length: Vec<usize>,
}

impl Simulation {
    fn initialize(n: usize, m: usize, arrival_mean: f64, service_mean: f64) -> Simulation {
        let mut servers = Vec::new();
        let mut lcgrand = LcGrand::new();
        for _ in 0..m {
            servers.push(Server::new());
        }
        let mut event_list = BinaryHeap::new();
        let arrival_time = lcgrand.expon(arrival_mean);
        event_list.push(Event::new(arrival_time, EventState::Arrival));

        Simulation {
            clock: 0.0,
            event_list,
            servers,
            lcgrand,
            arrival_mean,
            service_mean,
            queue: Vec::new(),
            client_count: n,
            served_clients: Vec::new(),
            queue_length: Vec::new(),
        }
    }

    fn clock_handler(&mut self) {
        let event = match self.event_list.pop() {
            Some(event) => event,
            None => panic!("Event list is empty"),
        };
        self.clock = event.time;
        match event.state {
            EventState::Arrival => self.arrival_handler(event),
            EventState::Departure(server_idx) => self.departure_handler(server_idx),
        }
    }

    fn arrival_handler(&mut self, event: Event) {
        self.queue.push(Client::new(event.time));
        self.queue_length.push(self.queue.len());
        self.client_count -= 1;

        for (i, server) in self.servers.iter_mut().enumerate() {
            match server.state {
                ServerState::Idle => {
                    self.service_handler(i);
                    break;
                }
                ServerState::Busy => {} // client must wait in queue
            }
        }

        if self.client_count == 0 {
            return;
        }

        // schedule next arrival
        let new_arrival_time = event.time + self.lcgrand.expon(self.arrival_mean);
        self.event_list
            .push(Event::new(new_arrival_time, EventState::Arrival));
    }

    fn departure_handler(&mut self, server: usize) {
        match self.servers[server].state {
            ServerState::Idle => panic!("Server is idle but departure event is scheduled"),
            ServerState::Busy => self.service_handler(server),
        }
    }

    fn service_handler(&mut self, server_idx: usize) {
        let mut client = match self.queue.pop() {
            Some(client) => client,
            None => {
                self.servers[server_idx].state = ServerState::Idle;
                return;
            }
        };
        client.service_start_time = self.clock;
        client.wait_time = client.service_start_time - client.arrival_time;
        client.service_time = self.lcgrand.expon(self.service_mean);
        client.departure_time = client.service_start_time + client.service_time;

        self.servers[server_idx].state = ServerState::Busy;
        self.servers[server_idx].active_time += client.service_time;
        self.event_list.push(Event::new(
            client.departure_time,
            EventState::Departure(server_idx),
        ));

        self.served_clients.push(client);
    }

    fn get_report(self) -> Report {
        let total_clients = self.served_clients.len();
        let total_servers = self.servers.len();
        let total_simulation_time = self.clock;
        let mut served_clients = self.served_clients;
        served_clients.sort_by(|a, b| a.arrival_time.partial_cmp(&b.arrival_time).unwrap());

        let mut mean_interarrival_time = 0.0;
        let mut mean_service_time = 0.0;
        let mut mean_wait_time = 0.0;

        for (client, next_client) in served_clients.iter().zip(served_clients.iter().skip(1)) {
            mean_interarrival_time += next_client.arrival_time - client.arrival_time;
            mean_service_time += client.service_time;
            mean_wait_time += client.wait_time;
        }
        mean_interarrival_time /= served_clients.len() as f64 - 1.0;
        mean_service_time /= served_clients.len() as f64;
        mean_wait_time /= served_clients.len() as f64;

        let mut servers_utilization = Vec::new();
        for server in self.servers {
            servers_utilization.push(server.active_time / self.clock);
        }

        let mut mean_queue_length = 0.0;
        let queue_length = self.queue_length.len();
        for queue_length in self.queue_length {
            mean_queue_length += queue_length as f64;
        }
        mean_queue_length /= queue_length as f64;

        Report {
            total_clients,
            total_servers,
            mean_wait_time,
            mean_service_time,
            mean_queue_length,
            total_simulation_time,
            mean_interarrival_time,
            servers_utilization,
        }
    }
}

pub fn calc_c_erlang(beta: f64, gamma: f64, m: f64, verbose: bool) -> f64 {
    fn factorial(n: i32) -> f64 {
        (1..=n).fold(1.0, |acc, x| acc * x as f64)
    }
    let lambda = 1.0 / beta;
    let mu = 1.0 / gamma;
    let rho = lambda / (m * mu);
    let lambda_over_mu = lambda / mu;
    let lambda_over_mu_to_m = lambda_over_mu.powf(m);
    let one_over_m_factorial = 1.0 / factorial(m as i32);
    let one_over_one_minus_rho = 1.0 / (1.0 - rho);
    let multiplier = lambda_over_mu_to_m * one_over_m_factorial * one_over_one_minus_rho;
    let mut sum = 0.0;
    for k in 1..m as i32 {
        sum += 1.0 / factorial(k) * lambda_over_mu.powi(k);
    }
    let p_0 = 1.0 / (1.0 + multiplier + sum);
    let c_erlang = p_0 * multiplier;

    if verbose {
        println!("---------C-Erlang---------");
        println!("{:<20}{:.3}", "β", beta);
        println!("{:<20}{:.3}", "γ", gamma);
        println!("{:<20}{:.3}", "m", m);
        println!("{:<20}{:.3}", "λ", lambda);
        println!("{:<20}{:.3}", "μ", mu);
        println!("{:<20}{:.3}", "ρ", rho);
        println!("{:<20}{:.3}", "λ/μ", lambda_over_mu);
        println!("{:<20}{:.3}", "λ/μ^m", lambda_over_mu_to_m);
        println!("{:<20}{:.3}", "1/m!", one_over_m_factorial);
        println!("{:<20}{:.3}", "1/(1-ρ)", one_over_one_minus_rho);
        println!("{:<20}{:.3}", "multiplier", multiplier);
        println!("{:<20}{:.3}", "sum", sum);
        println!("{:<20}{:.3}", "C-Erlang", c_erlang);
        println!("--------------------------");
    }

    c_erlang
}

fn main() {
    let args = Args::parse();
    let n = args.clients;
    let m = args.servers;
    let arrival_mean = args.arrival_mean;
    let service_mean = args.service_mean;
    let verbose = args.verbose;

    println!("Clients: {}", n);
    println!("Servers: {}", m);
    println!("Arrival Mean: {}", arrival_mean);
    println!("Service Mean: {}", service_mean);

    let mut simulation = Simulation::initialize(n, m, arrival_mean, service_mean);

    while simulation.client_count > 0 || !simulation.event_list.is_empty() {
        simulation.clock_handler();
    }

    let report = simulation.get_report();

    let teoric_c_erlang = calc_c_erlang(arrival_mean, service_mean, m as f64, verbose);
    let simulated_c_erlang = calc_c_erlang(
        report.mean_interarrival_time,
        report.mean_service_time,
        m as f64,
        verbose,
    );

    println!("Teoric C-Erlang: {:.4}", teoric_c_erlang);
    println!("Simulated C-Erlang: {:.4}", simulated_c_erlang);

    let error = (teoric_c_erlang - simulated_c_erlang).abs() / teoric_c_erlang;
    println!("Error: {:.3}%", error * 100.0);

    if verbose {
        println!();
        println!("Total Clients: {}", report.total_clients);
        println!("Total Servers: {}", report.total_servers);
        println!("Total Simulation Time: {:.3}", report.total_simulation_time);
        println!(
            "Mean Interarrival Time: {:.3}",
            report.mean_interarrival_time
        );
        println!("Mean Service Time: {:.3}", report.mean_service_time);
        println!("Mean Wait Time: {:.3}", report.mean_wait_time);
        println!("Mean Queue Length: {:.3}", report.mean_queue_length);

        for (i, server_utilization) in report.servers_utilization.iter().enumerate() {
            println!(
                "Server {} Utilization: {:.3}%",
                i + 1,
                server_utilization * 100.0
            );
        }
    }
}
