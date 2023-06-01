mod lcgrand;

use lcgrand::LcGrand;

#[derive(Clone, Debug, Eq, PartialEq)]
enum EventType {
    Free,
    Busy,
    Arrival,
}

#[derive(Clone, Debug)]
struct Event {
    event_type: EventType,
    time: f64,
    server: Option<i32>,
}

impl Event {
    fn new(server: Option<i32>) -> Self {
        match server {
            Some(_) => Event {
                event_type: EventType::Free,
                time: 1.0e+30,
                server,
            },
            None => Event {
                event_type: EventType::Arrival,
                time: 1.0e+30,
                server,
            },
        }
    }
}

impl Ord for Event {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.time.partial_cmp(&other.time).unwrap()
    }
}

impl PartialOrd for Event {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.time.partial_cmp(&other.time)
    }
}

impl PartialEq for Event {
    fn eq(&self, other: &Self) -> bool {
        self.time == other.time
    }
}

impl Eq for Event {}

struct Queue {
    lcgrand: LcGrand,
    events: Vec<Event>,
    waiting_clients: i32,
    median: f64,
}

impl Queue {
    fn new(num_servers: i32) -> Self {
        let mut events = Vec::new();
        events.push(Event::new(None));
        for i in 0..num_servers {
            events.push(Event::new(Some(i)));
        }
        Queue {
            lcgrand: LcGrand::new(),
            events,
            waiting_clients: 0,
            median: 1.0,
        }
    }

    fn handle_event(&mut self, time: f64) {
        self.events.sort();
        let event = self.events.first_mut().unwrap();
        match event.event_type {
            EventType::Arrival => {
                event.time = time + self.lcgrand.expon(self.median);
                if event.time == 1.0e+30 {
                    return;
                }

                let event = self
                    .events
                    .iter_mut()
                    .find(|event| event.event_type == EventType::Free);
                match event {
                    Some(event) => {
                        event.event_type = EventType::Busy;
                        event.time = time + self.lcgrand.expon(self.median);
                    }
                    None => {
                        self.waiting_clients += 1;
                    }
                }
            }
            EventType::Busy => {
                if self.waiting_clients > 0 {
                    self.waiting_clients -= 1;
                    event.time = time + self.lcgrand.expon(self.median);
                } else {
                    event.time = 1.0e+30;
                    event.event_type = EventType::Free;
                }
            }
            EventType::Free => {
                println!("Free event at time {:.2}", time);
                panic!("Free event should not be handled");
            }
        }
    }

    fn next_event(&self) -> Event {
        let mut events = self.events.clone();
        events.sort();
        events.first().unwrap().clone()
    }
}

fn main() {
    let mut queue = Queue::new(3);
    let mut time = 0.0;
    let mut clients = 100;

    while clients > 0 {
        queue.handle_event(time);
        let event = queue.next_event();
        println!("Event type: {:?}, time: {:.2}, server: {:?}", event.event_type, event.time, event.server);
        time = event.time;
        if event.event_type == EventType::Arrival {
            clients -= 1;
        }
    }
}
